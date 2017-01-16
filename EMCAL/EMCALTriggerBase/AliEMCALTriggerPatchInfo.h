#ifndef ALIEMCALTRIGGERPATCHINFO_H
#define ALIEMCALTRIGGERPATCHINFO_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TObject.h"

#include <TLorentzVector.h>
#include <TMath.h>
#include "AliEMCALTriggerBitConfig.h"
#include "AliEMCALTriggerConstants.h"


class AliEMCALGeometry;
class TArrayI;


/**
 * \class AliEMCALTriggerPatchInfo
 * \brief Main data structure storing all relevant information of EMCAL/DCAL trigger patches
 * \author Jiri Kral <>, University of Jyv&aumlskul&auml
 * \author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * \date Jun 26, 2013
 *
 * Emcal trigger patch information class
 * Can contain three types of information, distinguished by the various bits in the bit field:
 * -# online trigger information (no extra bits set)
 * -# offline recalculated trigger patches (bit 25, kSimpleOfflineTriggerBit set)
 * -# highest patch energy, also for events that did not fire the trigger (bits 22, 23 kRecalc... (using both online and offline info, use bit 25 to distinguish)
 */
class AliEMCALTriggerPatchInfo: public TObject {

 public:
  enum {
    kRecalcOffset = 16,
    kOfflineOffset = 24
  };

  enum CaloDetectorType_t {
    kEMCALdet = 0,
    kDCALPHOSdet = 1
  };

  /**
   * Default constructor
   */
  AliEMCALTriggerPatchInfo();

  /**
   * Copy constructor
   *
   * @param p Reference for the copy
   */
  AliEMCALTriggerPatchInfo(const AliEMCALTriggerPatchInfo &p);

  /**
   * Assignment operator
   *
   * @param p Reference for assignment
   * @return This object after assignment
   */
  AliEMCALTriggerPatchInfo &operator=(const AliEMCALTriggerPatchInfo &p);

  /**
   * Destructor
   */
  virtual ~AliEMCALTriggerPatchInfo();

  /**
   * Reset all fields to the default values using the standard constructor
   */
  void Reset();

  /**
   * Initialize patch
   * @param col0        Start column of the patch
   * @param row0        Start row of the patch
   * @param size        Size of the patch
   * @param adc         ADC signal of the patch
   * @param offlineAdc  Offline ADC signal of the patch
   * @param patchE      Energy of the patch (sum of cell amplitudes)
   * @param bitmask     Trigger bit mask of the patch
   * @param vertex      Primary vertex of the event
   * @param geom        Pointer to the EMCal geometry object
   */
  void Initialize(UChar_t col0, UChar_t row0, UChar_t size, UInt_t adc, UInt_t offlineAdc, Double_t patchE, UInt_t bitmask, const TVector3& vertex, const AliEMCALGeometry* geom);

  /**
   * Initialize patch
   * @param col0        Start column of the patch
   * @param row0        Start row of the patch
   * @param size        Size of the patch
   * @param adc         ADC signal of the patch
   * @param offlineAdc  Offline ADC signal of the patch
   * @param bitmask     Trigger bit mask of the patch
   * @param geom        Pointer to the EMCal geometry object
   */
  void Initialize(UChar_t col0, UChar_t row0, UChar_t size, UInt_t adc, UInt_t offlineAdc, UInt_t bitmask, const AliEMCALGeometry* geom);

  /**
   * Allocate a new AliEMCALTriggerPatchInfo object and initialize it
   * @param col0        Start column of the patch
   * @param row0        Start row of the patch
   * @param size        Size of the patch
   * @param adc         ADC signal of the patch
   * @param offlineAdc  Offline ADC signal of the patch
   * @param patchE      Energy of the patch (sum of cell amplitudes)
   * @param bitmask     Trigger bit mask of the patch
   * @param vertex      Primary vertex of the event
   * @param geom        Pointer to the EMCal geometry object
   * @return            Pointer to a new and initialized AliEMCALTriggerPatchInfo object (caller is responsible for releasing memory)
   */
  static AliEMCALTriggerPatchInfo* CreateAndInitialize(UChar_t col0, UChar_t row0, UChar_t size, UInt_t adc, UInt_t offlineAdc, Double_t patchE, UInt_t bitmask, const TVector3& vertex, const AliEMCALGeometry* geom);

  /**
   * Allocate a new AliEMCALTriggerPatchInfo object and initialize it
   * @param col0        Start column of the patch
   * @param row0        Start row of the patch
   * @param size        Size of the patch
   * @param adc         ADC signal of the patch
   * @param offlineAdc  Offline ADC signal of the patch
   * @param bitmask     Trigger bit mask of the patch
   * @param geom        Pointer to the EMCal geometry object
   * @return            Pointer to a new and initialized AliEMCALTriggerPatchInfo object (caller is responsible for releasing memory)
   */
  static AliEMCALTriggerPatchInfo* CreateAndInitialize(UChar_t col0, UChar_t row0, UChar_t size, UInt_t adc, UInt_t offlineAdc, UInt_t bitmask, const AliEMCALGeometry* geom);

  /**
   * Recalculate patch kinematic variables
   * @param patchE      Energy of the patch (sum of cell amplitudes)
   * @param vertex      Primary vertex of the event
   * @param geom        Pointer to the EMCal geometry object
   */
  void RecalculateKinematics(Double_t patchE, const TVector3& vertex, const AliEMCALGeometry* geom);

  /**
   * Access \f$ \phi \f$ angle of the geometric center of the trigger patch
   * @return \f$ \phi \f$ angle
   */
  Double_t GetPhiGeo() const { return GetPhiTransform(fCenterGeo.Phi()); }

  /**
   * Access \f$ \phi \f$ angle of the patch at the center of mass
   * @return \f$ \phi \f$ angle
   */
  Double_t GetPhiCM()  const { return GetPhiTransform(fCenterMass.Phi()); }

  /**
   * Get minimal \f$ \phi \f$ of the patch
   * @return \f$ \phi \f$ angle
   */
  Double_t GetPhiMin() const { return GetPhiTransform(fEdge1.Phi()); }

  /**
   * Get maximal \f$ \phi \f$ of the patch
   * @return \f$ \phi \f$ angle
   */
  Double_t GetPhiMax() const { return GetPhiTransform(fEdge2.Phi()); }

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
   * Return smeared patch energy. The smeared patch energy
   * is obtained from the offline energy smearing at each
   * FastOR. The smeared energy has to be provided by the
   * trigger maker or equivalent.
   * @return Smeared patch energy
   */
  Double_t GetSmearedEnergy() const { return fEnergySmeared; }

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
  Double_t GetADCAmpGeVRough() const { return (Double_t)fADCAmp * EMCALTrigger::kEMCL1ADCtoGeV; }

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

  /**
   * Return cell indices of the given patch in the cell array
   * @param geom EMCAL Geometry used in the run where the trigger patch was created from
   * @param cells Output array of cell indices corresponding to the given trigger patch
   */
  void     GetCellIndices( AliEMCALGeometry *geom, TArrayI *cells );
  
  /**
   * Get the transverse energy of the patch
   * @return Transverse energy of the patch
   */
  Double_t GetPatchET() const { return GetET(GetPatchE()); }

  /**
   * Get the online transverse energy of the patch
   * @return online transverse energy of the patch
   */
  Double_t GetPatchETfromADCAmp() const { return GetET(GetADCAmpGeVRough()); }

  /**
   * Get the starting row of the patch
   * @return Starting row of the patch
   */
  Int_t GetRowStart() const { return fRow0; }

  /**
   * Get the starting column of the patch
   * @return Starting column of the patch
   */
  Int_t GetColStart() const { return fCol0; }

  /**
   * Check whether patch is an EMCAL Level0 patch
   * @return True if patch is an EMCAL Level0 patch, false otherwise
   */
  Bool_t   IsLevel0() const { return TESTBIT(fTriggerBits, fOffSet + fTriggerBitConfig.GetLevel0Bit()); }

  /**
   * Check whether patch is an EMCAL Level1 jet patch passing the low threshold, found by the trigger electronics or the trigger simulation
   * @return True if patch matches the trigger condition, false otherwise
   */
  Bool_t   IsJetLow() const { return TESTBIT(fTriggerBits, fOffSet + fTriggerBitConfig.GetJetLowBit()); }

  /**
   * Check whether patch is an EMCAL Level1 jet patch passing the high threshold, found by the trigger electronics or the trigger simulation
   * @return True if patch matches the trigger condition, false otherwise
   */
  Bool_t   IsJetHigh() const { return TESTBIT(fTriggerBits, fOffSet + fTriggerBitConfig.GetJetHighBit()); }

  /**
   * Check whether patch is an EMCAL Level1 gamma patch passing the low threshold, found by the trigger electronics or the trigger simulation
   * @return True if patch matches the trigger condition, false otherwise
   */
  Bool_t   IsGammaLow() const { return TESTBIT(fTriggerBits, fOffSet + fTriggerBitConfig.GetGammaLowBit()); }

  /**
   * Check whether patch is an EMCAL Level1 gamma patch passing the high threshold, found by the trigger electronics or the trigger simulation
   * @return True if patch matches the trigger condition, false otherwise
   */
  Bool_t   IsGammaHigh() const { return TESTBIT(fTriggerBits, fOffSet + fTriggerBitConfig.GetGammaHighBit()); }

  /**
    * No background patches from hardware
    * @return Always false
    */
   Bool_t   IsBkg() const { return false; }

  /**
   * No main trigger any more
   * @return Always false
   */
  Bool_t   IsMainTrigger() const { return false; }

  /**
   * Check whether patch is an EMCAL Level1 jet patch passing the low threshold, found by the simple offline trigger
   * @return True if patch matches the trigger condition, false otherwise
   */
  Bool_t   IsLevel0Simple() const { return TESTBIT(fTriggerBits, kOfflineOffset + fTriggerBitConfig.GetLevel0Bit()); }

  /**
   * Check whether patch is an EMCAL Level1 jet patch passing the low threshold, found by the simple offline trigger
   * @return True if patch matches the trigger condition, false otherwise
   */
  Bool_t   IsJetLowSimple() const { return TESTBIT(fTriggerBits, kOfflineOffset + fTriggerBitConfig.GetJetLowBit()); }

  /**
   * Check whether patch is an EMCAL Level1 jet patch passing the high threshold, found by the simple offline trigger
   * @return True if patch matches the trigger condition, false otherwise
   */
  Bool_t   IsJetHighSimple() const { return TESTBIT(fTriggerBits, kOfflineOffset + fTriggerBitConfig.GetJetHighBit()); }

  /**
   * Check whether patch is an EMCAL Level1 gamma patch passing the low threshold, found by the simple offline trigger
   * @return True if patch matches the trigger condition, false otherwise
   */
  Bool_t   IsGammaLowSimple() const { return TESTBIT(fTriggerBits, kOfflineOffset + fTriggerBitConfig.GetGammaLowBit()); }

  /**
   * Check whether patch is an EMCAL Level1 gamma patch passing the high threshold, found by the simple offline trigger
   * @return True if patch matches the trigger condition, false otherwise
   */
  Bool_t   IsGammaHighSimple() const { return TESTBIT(fTriggerBits, kOfflineOffset + fTriggerBitConfig.GetGammaHighBit()); }

  /**
   * Check whether patch is an EMCAL background patch, found by the simple offline trigger
   * @return True if patch matches the trigger condition, false otherwise
   */
  Bool_t   IsBkgSimple() const { return TESTBIT(fTriggerBits, kOfflineOffset + fTriggerBitConfig.GetBkgBit()); }

  /**
   * Check whether patch is an EMCAL Level1 jet patch passing the low threshold, found by the simple offline trigger
   * @return True if patch matches the trigger condition, false otherwise
   */
  Bool_t   IsLevel0Recalc() const { return TESTBIT(fTriggerBits, kRecalcOffset + fTriggerBitConfig.GetLevel0Bit()); }

  /**
   * Check whether patch is an EMCAL Level1 jet patch passing the low threshold, found by the simple offline trigger
   * @return True if patch matches the trigger condition, false otherwise
   */
  Bool_t   IsJetLowRecalc() const { return TESTBIT(fTriggerBits, kRecalcOffset + fTriggerBitConfig.GetJetLowBit()); }

  /**
   * Check whether patch is an EMCAL Level1 jet patch passing the high threshold, found by the simple offline trigger
   * @return True if patch matches the trigger condition, false otherwise
   */
  Bool_t   IsJetHighRecalc() const { return TESTBIT(fTriggerBits, kRecalcOffset + fTriggerBitConfig.GetJetHighBit()); }

  /**
   * Check whether patch is an EMCAL Level1 gamma patch passing the low threshold, found by the simple offline trigger
   * @return True if patch matches the trigger condition, false otherwise
   */
  Bool_t   IsGammaLowRecalc() const { return TESTBIT(fTriggerBits, kRecalcOffset + fTriggerBitConfig.GetGammaLowBit()); }

  /**
   * Check whether patch is an EMCAL Level1 gamma patch passing the high threshold, found by the simple offline trigger
   * @return True if patch matches the trigger condition, false otherwise
   */
  Bool_t   IsGammaHighRecalc() const { return TESTBIT(fTriggerBits, kRecalcOffset + fTriggerBitConfig.GetGammaHighBit()); }

  /**
   * Check whether patch is an EMCAL background patch, found by the simple offline trigger
   * @return True if patch matches the trigger condition, false otherwise
   */
  Bool_t   IsBkgRecalc() const { return TESTBIT(fTriggerBits, kRecalcOffset + fTriggerBitConfig.GetBkgBit()); }

  /**
   * Now main trigger any more in the new definition
   * @return Always false
   */
  Bool_t   IsMainTriggerSimple() const { return kFALSE; }

  /**
   * Check whether patch is found by the simple offline trigger (on offline amplitudes)
   * @return True if the patch is found by the simple offline trigger, false otherwise
   */
  Bool_t   IsOfflineSimple() const { return IsJetLowSimple() || IsJetHighSimple() || IsGammaLowSimple() || IsGammaHighSimple(); }

  /**
   * Check whether patch is found by the simple offline trigger (on offline amplitudes)
   * @return True if the patch is found by the simple offline trigger, false otherwise
   */
  Bool_t   IsRecalc() const { return IsJetLowRecalc() || IsJetHighRecalc() || IsGammaLowRecalc() || IsGammaHighRecalc(); }

  /**
   * Check whether patch is found by the simple offline trigger (on offline amplitudes)
   * @return True if the patch is found by the simple offline trigger, false otherwise
   */
  Bool_t   IsOnline() const { return IsJetLow() || IsJetHigh() || IsGammaLow() || IsGammaHigh(); }

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
  Bool_t   IsRecalcJet() const { return IsJetLowRecalc() || IsJetHighRecalc(); }

  /**
   * Check if the patch is a recalculated gamma patch
   * @return True if the patch is a recalculated gamma patch, false otherwise
   */
  Bool_t   IsRecalcGamma() const { return IsGammaLowRecalc() || IsGammaHighRecalc(); }

  /**
   * Check whether patch is in the EMCal
   * @return True if the patch is in the EMCal
   */
  Bool_t         IsEMCal()         const { return (fDetectorType == kEMCALdet)   ; }

  /**
   * Check whether patch is in the EMCal
   * @return True if the patch is in the EMCal
   */
  Bool_t         IsDCalPHOS()      const { return (fDetectorType == kDCALPHOSdet); }

  /**
   * @return Detector in which the patch is located (EMCal or DCal/PHOS)
   */
  CaloDetectorType_t GetDetectorType() const { return fDetectorType                  ; }

  /**
   * Check whether a trigger bit is set
   * @param bitnumber Bit number to be tested
   * @return True if the bit is set
   */
  Bool_t   TestTriggerBit(UInt_t bitnumber) const { return TESTBIT(fTriggerBits, bitnumber); }

  /**
   * Returns the patch size
   * @return patch size
   */
  UChar_t  GetPatchSize() const { return fPatchSize; }

  /**
   * Set the starting row
   * @param row0 Starting row of the patch
   */
  void SetRowStart(int row0) { fRow0 = row0; }

  /**
   * Set the starting column
   * @param col0 Starting column of the patch
   */
  void SetCol0(int col0) { fCol0 = col0; }

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
   * Set the smeared energy. Smeared energy is defined as energy obtained
   * from the offline energy by randomizing the energy with the appropriate
   * distribution (gauss with proper values for mean and width).
   * @param e
   */
  void SetSmearedEnergy(Double_t e) { fEnergySmeared = e; }

  /**
   * Set Indices in x and y of the edge cell
   * @param x Cell index in x-direction
   * @param y Cell index in y-direction
   */
  void SetEdgeCell( Int_t x, Int_t y ) { fEdgeCell[0] = x; fEdgeCell[1] = y; }

  /**
   * No simple offline bit any more as patches contain combined information from
   * online, offline and recalc trigger
   */
  void SetOfflineSimple() { }

  /**
   * Define Lorentz vector of the given trigger patch
   * @param lv Lorentz vector to be defined
   * @param v Patch vector position
   * @param e Patch energy
   */
  void SetLorentzVector( TLorentzVector &lv, const TVector3 &v, Double_t e );

  /**
   * Set the trigger bits
   * @param i Trigger bits of the patch
   */
  void SetTriggerBits( Int_t i ) { fTriggerBits = i; }

  /**
   * Set detector in which the patch is located (EMCal or DCal/PHOS)
   * @param t Detector type
   */
  void SetDetectorType(CaloDetectorType_t t)  { fDetectorType = t; }

  /**
   * Set the MC trigger bit offset
   * @param i MC trigger bit offset
   */
  void SetOffSet(Int_t i)        { fOffSet      = i; }

  /**
   * Set the trigger bit configuration
   * @param ref Trigger bit configuration used to create the patch
   */
  void SetTriggerBitConfig(const AliEMCALTriggerBitConfig * ref) { fTriggerBitConfig.Initialise(*ref); }

  /**
   * Get the trigger bit configuration used to create the trigger patch
   * @return Trigger bit configuration of the patch
   */
  const AliEMCALTriggerBitConfig *GetTriggerBitConfig() const { return &fTriggerBitConfig; }


 protected:
  Double_t GetPhiTransform(Double_t phiin) const;
  Double_t GetET(Double_t energy) const;

  TLorentzVector    fCenterGeo;                     ///< geometrical center
  TLorentzVector    fCenterMass;                    ///< CM
  TLorentzVector    fEdge1;                         ///< max eta/ min phi edge
  TLorentzVector    fEdge2;                         ///< min eta/ max phi edge
  Int_t             fADCAmp;                        ///< online (trigger) ADC amplitude
  Int_t             fADCOfflineAmp;                 ///< offline (FEE) ADC amplitude
  Double_t          fEnergySmeared;                 ///< Smeared patch energy
  Int_t             fTriggerBits;                   ///< trigger bit mask, see definitions in AliEMCALTriggerType and TriggerMakerBits_t (above)
  Int_t             fEdgeCell[2];                   ///< cell "bottom lower" edge (min phi, max eta)
  Int_t             fOffSet;                        ///< offset of bit (different in data and MC)
  Int_t             fCol0;                          ///< Start column
  Int_t             fRow0;                          ///< Start row
  UChar_t           fPatchSize;                     ///< Trigger patch size
  CaloDetectorType_t    fDetectorType;                  ///< Detector type (EMCal or DCal/PHOS)
  AliEMCALTriggerBitConfig fTriggerBitConfig;     ///< Trigger bit configuration

  /// \cond CLASSIMP
  ClassDef(AliEMCALTriggerPatchInfo, 6) // Emcal particle class
  /// \endcond
};
#endif
