/**
 * \file AliEmcalTriggerPatchInfo.h
 * \brief Class to make array of trigger patch objects in AOD/ESD events.
 *
 * Emcal trigger patch information class
 * Can contain three types of information, distinguished by the various bits in the bit field:
 * -# online trigger information (no extra bits set)
 * -# offline recalculated trigger patches (bit 25, kSimpleOfflineTriggerBit set)
 * -# highest patch energy, also for events that did not fire the trigger (bits 22, 23 kRecalc... (using both online and offline info, use bit 25 to distinguish)
 *
 * \author Jiri Kral <>, University of Jyv&aumlskul&auml
 * \date Jun 26, 2013
 */
#ifndef ALIEMCALTRIGGERPATCHINFO_H
#define ALIEMCALTRIGGERPATCHINFO_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TObject.h"

#include <TLorentzVector.h>
#include <TMath.h>
#include "AliEmcalTriggerBitConfig.h"
#include "AliEmcalTriggerSetupInfo.h"

class AliEMCALGeometry;
class TArrayI;

/**
 * \class AliEmcalTriggerPatchInfo
 * \brief Main data structure storing all relevant information of EMCAL/DCAL trigger patches
 *
 * Emcal trigger patch information class
 * Can contain three types of information, distinguished by the various bits in the bit field:
 * -# online trigger information (no extra bits set)
 * -# offline recalculated trigger patches (bit 25, kSimpleOfflineTriggerBit set)
 * -# highest patch energy, also for events that did not fire the trigger (bits 22, 23 kRecalc... (using both online and offline info, use bit 25 to distinguish)
 */
class AliEmcalTriggerPatchInfo: public TObject {
 public:
  AliEmcalTriggerPatchInfo();
  AliEmcalTriggerPatchInfo(const AliEmcalTriggerPatchInfo &p); 
  AliEmcalTriggerPatchInfo &operator=(const AliEmcalTriggerPatchInfo &p);
  virtual ~AliEmcalTriggerPatchInfo();

  /**
   * \enum Bit definition for special trigger bits
   */
  enum TriggerMakerBits_t {
    kRecalcJetBitNum = 22,   ///< Trigger bit for recalculated jet patches
    kRecalcGammaBitNum = 23, ///< Trigger bit for recalculated gamma patches
    kMainTriggerBitNum = 24, ///< Trigger bit indicating the main (highest energy) trigger patch of a given type per event
    kSimpleOfflineBitNum = 25///< Trigger bit indicating that the patch was created by the offline trigger algorithm
  };

  /**
   * Access \f$ \phi \f$ angle of the geometric center of the trigger patch
   * @return \f$ \phi \f$ angle
   */
  Double_t GetPhiGeo() const { return fCenterGeo.Phi(); }
  /**
   * Access \f$ \phi \f$ angle of the patch at the center of mass
   * @return \f$ \phi \f$ angle
   */
  Double_t GetPhiCM()  const { return fCenterMass.Phi(); }
  /**
   * Get minimal \f$ \phi \f$ of the patch
   * @return \f$ \phi \f$ angle
   */
  Double_t GetPhiMin() const { return fEdge1.Phi(); }
  /**
   * Get maximal \f$ \phi \f$ of the patch
   * @return \f$ \phi \f$ angle
   */
  Double_t GetPhiMax() const { return fEdge2.Phi(); }
  /**
   * Get \f$ \eta \f$ of the patch at the geometrical center
   * @return Patch \f$ \eta \f$
   */
  Double_t GetEtaGeo() const { return fCenterGeo.Eta(); }
  /**
   * Get \f$ \eta \f$ of the patch at the center of mass
   * @return Patch \f$ \eta \f$
   */
  Double_t GetEtaCM()  const { return fCenterMass.Eta(); }
  /**
   * Get minimum \f$ \eta \f$ of the patch
   * @return Patch \f$ \eta \f$
   */
  Double_t GetEtaMin() const { return fEdge2.Eta(); }
  /**
   * Get maximum \f$ \eta \f$ of the patch
   * @return Patch \f$ \eta \f$
   */
  Double_t GetEtaMax() const { return fEdge1.Eta(); }
  /**
   * Get the patch energy
   * @return Patch energy
   */
  Double_t GetPatchE() const { return fCenterGeo.E(); }
  /**
   * Get patch online ADC amplitude
   * @return Online patch amplitude
   */
  Int_t    GetADCAmp() const { return fADCAmp; }
  /**
   * Get patch offline ADC amplitude (obtained from calibrated cell energies converted into ADC signals)
   * @return Offline patch amplitude
   */
  Int_t    GetADCOfflineAmp() const { return fADCOfflineAmp; }
  /**
   * Get patch energy estimated from offline ADC amplitude converted into energya
   * @return Patch energy estimate
   */
  Double_t GetADCAmpGeVRough() const { return (Double_t)fADCAmp * kEMCL1ADCtoGeV; }
  /**
   * Get the trigger bits of the classes which fired the patch
   * @return Selected trigger bits
   */
  Int_t    GetTriggerBits() const { return fTriggerBits; }
  /**
   * Get X position of the edge cell
   * @return Cell x-position
   */
  Int_t    GetEdgeCellX() const { return fEdgeCell[0]; }
  /**
   * Get Y position of the edge cell
   * @return Cell y-position
   */
  Int_t    GetEdgeCellY() const { return fEdgeCell[1]; }
  void     GetCellIndices( AliEMCALGeometry *geom, TArrayI *cells );
  
  /**
   * Check whether patch is an EMCAL Level0 patch
   * @return True if patch is an EMCAL Level0 patch, false otherwise
   */
  Bool_t   IsLevel0() const { return (Bool_t)((fTriggerBits >> (fOffSet + fTriggerBitConfig.GetLevel0Bit()))&(!(fTriggerBits >> kSimpleOfflineBitNum))&1); }
  /**
   * Check whether patch is an EMCAL Level1 jet patch passing the low threshold, found by the trigger electronics or the trigger simulation
   * @return True if patch matches the trigger condition, false otherwise
   */
  Bool_t   IsJetLow() const { return (Bool_t)((fTriggerBits >> (fOffSet + fTriggerBitConfig.GetJetLowBit()))&(!(fTriggerBits >> kSimpleOfflineBitNum))&1); }
  /**
   * Check whether patch is an EMCAL Level1 jet patch passing the high threshold, found by the trigger electronics or the trigger simulation
   * @return True if patch matches the trigger condition, false otherwise
   */
  Bool_t   IsJetHigh() const { return (Bool_t)((fTriggerBits >> (fOffSet + fTriggerBitConfig.GetJetHighBit()))&(!(fTriggerBits >> kSimpleOfflineBitNum))&1); }
  /**
   * Check whether patch is an EMCAL Level1 gamma patch passing the low threshold, found by the trigger electronics or the trigger simulation
   * @return True if patch matches the trigger condition, false otherwise
   */
  Bool_t   IsGammaLow() const { return (Bool_t)((fTriggerBits >> (fOffSet + fTriggerBitConfig.GetGammaLowBit()))&(!(fTriggerBits >> kSimpleOfflineBitNum))&1); }
  /**
   * Check whether patch is an EMCAL Level1 gamma patch passing the high threshold, found by the trigger electronics or the trigger simulation
   * @return True if patch matches the trigger condition, false otherwise
   */
  Bool_t   IsGammaHigh() const { return (Bool_t)((fTriggerBits >> (fOffSet + fTriggerBitConfig.GetGammaHighBit()))&(!(fTriggerBits >> kSimpleOfflineBitNum))&1); }
  /**
   * Check whether patch is the main EMCAL trigger patch of a given trigger type, found by the trigger electronics or the trigger simulation
   * @return True if patch is the main trigger patch, false otherwise
   */
  Bool_t   IsMainTrigger() const { return (Bool_t)((fTriggerBits >> kMainTriggerBitNum)&(!(fTriggerBits >> kSimpleOfflineBitNum))&1); }
  /**
   * Check whether patch is an EMCAL Level1 jet patch passing the low threshold, found by the simple offline trigger
   * @return True if patch matches the trigger condition, false otherwise
   */
  Bool_t   IsJetLowSimple() const { return (Bool_t)((fTriggerBits >> (fOffSet + fTriggerBitConfig.GetJetLowBit()))&(fTriggerBits >> kSimpleOfflineBitNum)&1); }
  /**
   * Check whether patch is an EMCAL Level1 jet patch passing the high threshold, found by the simple offline trigger
   * @return True if patch matches the trigger condition, false otherwise
   */
  Bool_t   IsJetHighSimple() const { return (Bool_t)((fTriggerBits >> (fOffSet + fTriggerBitConfig.GetJetHighBit()))&(fTriggerBits >> kSimpleOfflineBitNum)&1); }
  /**
   * Check whether patch is an EMCAL Level1 gamma patch passing the low threshold, found by the simple offline trigger
   * @return True if patch matches the trigger condition, false otherwise
   */
  Bool_t   IsGammaLowSimple() const { return (Bool_t)((fTriggerBits >> (fOffSet + fTriggerBitConfig.GetGammaLowBit()))&(fTriggerBits >> kSimpleOfflineBitNum)&1); }
  /**
   * Check whether patch is an EMCAL Level1 gamma patch passing the high threshold, found by the simple offline trigger
   * @return True if patch matches the trigger condition, false otherwise
   */
  Bool_t   IsGammaHighSimple() const { return (Bool_t)((fTriggerBits >> (fOffSet + fTriggerBitConfig.GetGammaHighBit()))&(fTriggerBits >> kSimpleOfflineBitNum)&1); }
  /**
   * Check whether patch is the main EMCAL trigger patch of a given trigger type, found by the simple offline trigger
   * @return True if patch is the main trigger patch, false otherwise
   */
  Bool_t   IsMainTriggerSimple() const { return (Bool_t)((fTriggerBits >> kMainTriggerBitNum)&(fTriggerBits >> kSimpleOfflineBitNum)&1); }
  /**
   * Check whether patch is found by the simple offline trigger (on offline amplitudes)
   * @return True if the patch is found by the simple offline trigger, false otherwise
   */
  Bool_t   IsOfflineSimple() const { return (Bool_t)((fTriggerBits >> kSimpleOfflineBitNum)&1); }

  /**
   * Access to Lorentz Vector of the centre-of-mass of the trigger patch
   * @return Lorentz Vector of the centre-of-mass of the trigger patch
   */
  const TLorentzVector &GetLorentzVectorCM() const { return fCenterMass; }
  /**
   * Access to Lorentz Vector of the geometrical centre of the trigger patch
   * @return Lorentz Vector of the geometrical centre of the trigger patch
   */
  const TLorentzVector &GetLorentzVectorCenterGeo() const { return fCenterGeo; }
  /**
   * Access to Lorentz Vector of the lower left edge of the trigger patch
   * @return Lorentz Vector of the lower left edge of the trigger patch
   */
  const TLorentzVector &GetLorentzVectorEdge1() const { return fEdge1; }
  /**
   * Access to Lorentz Vector of the upper right edge of the trigger patch
   * @return Lorentz Vector of the upper right edge of the trigger patch
   */
  const TLorentzVector &GetLorentzVectorEdge2() const { return fEdge2; }

  // Recalculated max patches
  /**
   * Check if the patch is a recalculated jet patch
   * @return True if the patch is a recalculated jet patch, false otherwise
   */
  Bool_t   IsRecalcJet() const { return (Bool_t) ((fTriggerBits >> kRecalcJetBitNum)&1); }
  /**
   * Check if the patch is a recalculated gamma patch
   * @return True if the patch is a recalculated gamma patch, false otherwise
   */
  Bool_t   IsRecalcGamma() const { return (Bool_t) ((fTriggerBits >> kRecalcGammaBitNum)&1); }

 /**
  * Set the geometric center position of the patch
  * @param v Position 3-vector
  * @param e Patch energy
  */
  void SetCenterGeo( const TVector3 &v, Double_t e ) { SetLorentzVector( fCenterGeo, v, e ); }
  /**
   * Set the geometric center position of the patch
   * @param v Position Lorentz vector
   */
  void SetCenterGeo( const TLorentzVector &v ) { fCenterGeo = v; }
  /**
   * Set the center-of-mass position of the trigger patch
   * @param v Position Lorentz vector
   */
  void SetCenterMass( const TLorentzVector &v ) { fCenterMass = v; }
  /**
   * Set the center-of-mass position of the trigger patch
   * @param v Position 3-vector
   * @param e Patch energy
   */
  void SetCenterMass( const TVector3 &v, Double_t e ) { SetLorentzVector( fCenterMass, v, e ); }
  /**
   * Set lower edge position of the trigger patch
   * @param v Position Lorentz vector
   */
  void SetEdge1( const TLorentzVector &v ) { fEdge1 = v; }
  /**
   * Set lower edge position of the trigger patch
   * @param v Position 3-vector
   * @param e Patch energy
   */
  void SetEdge1( const TVector3 &v, Double_t e ) { SetLorentzVector( fEdge1, v, e ); }
  /**
   * Set upper edge position of the trigger patch
   * @param v Lorentz-vector of the upper edge position of the trigger patch
   */
  void SetEdge2( const TLorentzVector &v ) { fEdge2 = v; }
  /**
   * Set upper edge position of the trigger patch
   * @param v Position 3-vector
   * @param e Patch Energy
   */
  void SetEdge2( const TVector3 &v, Double_t e ) { SetLorentzVector( fEdge2, v, e ); }
  /**
   * Set online ADC amplitude
   * @param a Online ADC amplitude
   */
  void SetADCAmp( Int_t a ) { fADCAmp = a; }
  /**
   * Set offline ADC amplitude (derived from cell energies converted to ADC amplitude)
   * @param a Offline ADC amplitude
   */
  void SetADCOfflineAmp( Int_t a ) { fADCOfflineAmp = a; }
  /**
   * Set Indices in x and y of the edge cell
   * @param x Cell index in x-direction
   * @param y Cell index in y-direction
   */
  void SetEdgeCell( Int_t x, Int_t y ) { fEdgeCell[0] = x; fEdgeCell[1] = y; }
  /**
   * Mark patch as created by the simple offline trigger
   */
  void SetOfflineSimple() { fTriggerBits |= 1 << kSimpleOfflineBitNum; }

  void SetLorentzVector( TLorentzVector &lv, const TVector3 &v, Double_t e );

  /**
   * Set the trigger bits
   * @param i Trigger bits of the patch
   */
  void SetTriggerBits( Int_t i ) { fTriggerBits = i; }

  /**
   * Set the MC trigger bit offset
   * @param i MC trigger bit offset
   */
  void SetOffSet(Int_t i)        { fOffSet      = i; }

  /**
   * Set the trigger bit configuration
   * @param ref Trigger bit configuration used to create the patch
   */
  void SetTriggerBitConfig(const AliEmcalTriggerBitConfig * ref) { fTriggerBitConfig.Initialise(*ref); }

  /**
   * Get the trigger bit configuration used to create the trigger patch
   * @return Trigger bit configuration of the patch
   */
  const AliEmcalTriggerBitConfig *GetTriggerBitConfig() const { return &fTriggerBitConfig; }


 protected:
  //TLorentzVector   &GetLorentzVector(const Double_t *vertex = 0)  const;

  TLorentzVector    fCenterGeo;                     ///< geometrical center
  TLorentzVector    fCenterMass;                    ///< CM
  TLorentzVector    fEdge1;                         ///< max eta/ min phi edge
  TLorentzVector    fEdge2;                         ///< min eta/ max phi edge
  Int_t             fADCAmp;                        ///< online (trigger) ADC amplitude
  Int_t             fADCOfflineAmp;                 ///< offline (FEE) ADC amplitude
  Int_t             fTriggerBits;                   ///< trigger bit mask, see definitions in AliEmcalTriggerType and TriggerMakerBits_t (above)
  Int_t             fEdgeCell[2];                   ///< cell "bottom lower" edge (min phi, max eta)
  Int_t             fOffSet;                        ///< offset of bit (different in data and MC)
  AliEmcalTriggerBitConfig   fTriggerBitConfig;     ///< Trigger bit configuration

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerPatchInfo, 6) // Emcal particle class
  /// \endcond
};
#endif
