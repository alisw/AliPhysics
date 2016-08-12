#ifndef ALIFMDRECONSTRUCTOR_H
#define ALIFMDRECONSTRUCTOR_H
//
//  Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
//  reserved. 
//
//  See cxx source for full Copyright notice                               
//
//  AliFMDReconstructor.h 
//  Task Class for making TreeR for FMD                        
//
//-- Authors: Evgeny Karpechev (INR) and Alla Maevskaia (INR)
//   Latest changes by Christian Holm Christensen <cholm@nbi.dk>
/* $Id$ */
/** @file    AliFMDReconstructor.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:47:09 2006
    @brief   FMD reconstruction 
*/

//____________________________________________________________________
// Header guards in the header files speeds up the compilation
// considerably.  Please leave them in. 
#ifndef ALIRECONSTRUCTOR_H
# include <AliReconstructor.h>
#endif
#include "AliLog.h"
#include <AliFMDBoolMap.h>

//____________________________________________________________________
class TTree;
class TClonesArray;
class AliFMDDigit;
class AliRawReader;
class AliFMDRawReader;
class AliESDEvent;
class AliESDFMD;
class AliFMDRecoParam;
class TH1;


/** 
 * @defgroup FMD_rec Reconstruction 
 *
 * Classes used for reconstruction of FMD data 
 *
 * @ingroup FMD 
 */
//____________________________________________________________________
/** 
 * @brief This is a class that reconstructs AliFMDRecPoint objects
 *        from of Digits.  
 *
 * This class reads either digits from a TClonesArray or raw data
 * from a DDL file (or similar), and applies calibrations to get
 * psuedo-inclusive multiplicities per strip.
 * 
 * @ingroup FMD_rec
 */
class AliFMDReconstructor: public AliReconstructor 
{
public:
  /** 
   * CTOR 
   */
  AliFMDReconstructor();
  /** 
   * DTOR 
   */
  virtual ~AliFMDReconstructor();

  /** 
   * Initialize the reconstructor.  Here, we initialize the geometry
   * manager, and finds the local to global transformations from the
   * geometry.   The calibration parameter manager is also
   * initialized (meaning that the calibration parameters is read
   * from CDB).
   */
  virtual void   Init();
  /** 
   * Flag that we can convert raw data into digits. 
   *
   * @return always @c true 
   */
  virtual Bool_t HasDigitConversion() const { return kTRUE; }
  /** 
   * Convert raw data read from the AliRawReader @a reader into
   * digits.  This is done using AliFMDRawReader and
   * AliFMDAltroReader.  The digits are put in the passed TTree @a
   * digitsTree. 
   *
   * @note This is the first part of the reconstruction as done by the
   * offical steering class AliReconstruction.
   *
   * @param reader Raw reader.  @param digitsTree Tree to store read
   * digits in.
   */
  virtual void   ConvertDigits(AliRawReader* reader, TTree* digitsTree) const;
  /** 
   * Reconstruct one event from the digits passed in @a digitsTree.
   * The member function creates AliFMDRecPoint objects and stores
   * them on the output tree @a clusterTree.  An FMD ESD object is
   * created in parallel. 
   *
   * @note This is the second part of the reconstruction as done by
   * the offical steering class AliReconstruction.
   *
   * @param digitsTree  Tree holding the digits of this event
   * @param clusterTree Tree to store AliFMDRecPoint objects in. 
   */
  virtual void   Reconstruct(TTree* digitsTree, TTree* clusterTree) const;
  /** 
   * Not used 
   * @todo Implement this, such that we'll reconstruct directly from
   *       the read ADC values rather than going via an intermedant
   *       TClonesArray of AliFMDDigits
   */
  virtual void   Reconstruct(AliRawReader *, TTree*) const;
  /** 
   * Not used.
   *
   * @todo This is called by the above same member function but with a
   * pointer to a AliRawReader object and a pointer to a TTree object.
   *
   * @param reader Reader object 
   */
  virtual void Reconstruct(AliFMDRawReader& reader) const;
  /** 
   * Put in the ESD data, the FMD ESD data.  The object created by
   * the Reconstruct member function is copied to the ESD object. 
   *
   * @note This is the third part of the reconstruction as done by
   * the offical steering class AliReconstruction.
   *
   * @param digitsTree   Tree of digits for this event - not used
   * @param clusterTree  Tree of reconstructed points for this event -
   *        not used.
   * @param esd ESD object to store data in. 
   */
  virtual void   FillESD(TTree* digitsTree, TTree* clusterTree, 
			 AliESDEvent* esd) const;
  /** 
   * Forwards to above member function 
   */
  virtual void   FillESD(AliRawReader*, TTree* clusterTree, 
			 AliESDEvent* esd) const;
  /** 
   * Return the filled FMD ESD object
   * 
   * @return FMD ESD object
   */
  AliESDFMD* GetESDObject() const { return fESDObj; }
  /** 
   * Create SDigits from raw data
   * 
   * @param reader  The raw reader
   * @param sdigits Array to fill with AliFMDSDigit objects. 
   */  
  virtual void Digitize(AliRawReader* reader, 
			TClonesArray* sdigits) const;
  
  /** 
   * Not used 
   */
  virtual void   SetESD(AliESDEvent* esd) { fESD = esd; }
  /** 
   * Set the noise factor 
   *
   * @param f Factor to use 
   */
  virtual void SetNoiseFactor(Float_t f=3) { fNoiseFactor = f; }
  /** 
   * Set whether we should do angle correction or nor 
   *
   * @param use If true, do angle correction 
   */
  virtual void SetAngleCorrect(Bool_t use=kTRUE) { fAngleCorrect = use; }
  /** 
   * Set whether we want to do diagnostics.  If this is enabled, a
   * file named @c FMD.Diag.root will be made.  It contains a set of
   * histograms for each event, filed in separate directories in the
   * file.  The histograms are 
   * @verbatim 
   * diagStep1   Read ADC vs. Noise surpressed ADC 
   * diagStep2   Noise surpressed ADC vs. calculated Energy dep.
   * diagStep3   Energy deposition vs. angle corrected Energy dep.
   * diagStep4   Energy deposition vs. calculated multiplicity
   * diagAll     Read ADC vs. calculated multiplicity
   * @endverbatim 
   *
   * @param use If true, make the diagnostics file 
   */
  void SetDiagnose(Bool_t use=kTRUE) { fDiagnostics = use; }
  /** 
   * Process AliFMDDigit objects in @a digits.  For each digit, find
   * the psuedo-rapidity @f$ \eta@f$, azimuthal angle @f$ \varphi@f$,
   * energy deposited @f$ E@f$, and psuedo-inclusive multiplicity @f$
   * M@f$.
   * 
   * @param digits  Array of digits. 
   * @param rawRead Raw reader used 
   */
  virtual void ProcessDigits(TClonesArray* digits,
			     const AliFMDRawReader& rawRead) const;
    
protected:
  /** 
   * Copy CTOR 
   *
   * @param other Object to copy from. 
   */
  AliFMDReconstructor(const AliFMDReconstructor&); //Not implemented
  /** 
   * Assignment operator 
   *
   * @param other Object to assign from
   *
   * @return reference to this object 
   */
  AliFMDReconstructor& operator=(const AliFMDReconstructor&); //Not implemented
  /** 
   * Run some checks before reconstruction, clear internal arrays, etc. 
   * 
   * @return true on success 
   */
  Bool_t PreReconstruct() const;
  /** 
   * Try to get the vertex from either ESD or generator header.  Sets
   * @c fCurrentVertex to the found Z posistion of the vertex (if 
   * found), and sets the flag @c fVertexType accordingly 
   *
   * @param esd ESD structure to get Vz from
   */
  virtual void GetVertex(AliESDEvent* esd) const;
  /** 
   * Set-up reconstructor to use values from reconstruction
   * parameters, if present, for this event.   If the argument @a set
   * is @c false, then restore preset values. 
   * 
   * @param set 
   */  
  virtual void UseRecoParam(Bool_t set=kTRUE) const;
  /** 
   * Process AliFMDDigit objects in @a digits.  For each digit, find
   * the psuedo-rapidity @f$ \eta@f$, azimuthal angle @f$ \varphi@f$,
   * energy deposited @f$ E@f$, and psuedo-inclusive multiplicity @f$
   * M@f$.
   * 
   * @param digits Array of digits. 
   */
  virtual void ProcessDigits(TClonesArray* digits) const;
  /** 
   * Process a single digit 
   * 
   * @param digit Digiti to process
   */ 
  virtual void ProcessDigit(AliFMDDigit* digit) const;
  /** 
   * Process the signal from a single strip. 
   * 
   * @param det Detector number 
   * @param rng Ring identifier 
   * @param sec Sector number
   * @param str Strip number 
   * @param adc Number of ADC counts for this strip
   */  
  virtual void ProcessSignal(UShort_t det, 
			     Char_t   rng, 
			     UShort_t sec, 
			     UShort_t str, 
			     Short_t  adc) const;
  /** 
   * Process the signal from a single strip. 
   * 
   * @param sdigits Array to fill
   * @param det     Detector number 
   * @param rng     Ring identifier 
   * @param sec     Sector number
   * @param str     Strip number 
   * @param sam     Sample number 
   * @param adc     Number of ADC counts for this strip
   */  
  virtual void DigitizeSignal(TClonesArray* sdigits, 
			      UShort_t      det, 
			      Char_t         rng, 
			      UShort_t       sec, 
			      UShort_t       str, 
			      UShort_t       sam,
			      Short_t        adc) const;
  /** 
   * Subtract the pedestal off the ADC counts. 
   * 
   * @param det           Detector number
   * @param rng           Ring identifier
   * @param sec           Sector number
   * @param str           Strip number
   * @param adc           ADC counts
   * @param noiseFactor   If pedestal substracted pedestal is less then
   *        this times the noise, then consider this to be 0. 
   * @param zsEnabled     Whether zero-suppression is on.
   * @param zsNoiseFactor Noise factor used in on-line pedestal
   *        subtraction. 
   * 
   * @return The pedestal subtracted ADC counts (possibly 0), or @c
   *         USHRT_MAX in case of problems.
   */  
  virtual UShort_t SubtractPedestal(UShort_t det, 
				    Char_t   rng, 
				    UShort_t sec, 
				    UShort_t str, 
				    UShort_t adc, 
				    Float_t  noiseFactor,
				    Bool_t   zsEnabled, 
				    UShort_t zsNoiseFactor) const;
  /** 
   * Substract pedestals from raw ADC in @a digit
   * 
   * @param det	  Detector number  
   * @param rng   Ring identifier 
   * @param sec   Sector number
   * @param str   Strip number 
   * @param adc   Number of ADC counts
   *
   * @return Pedestal subtracted ADC count. 
   */
  virtual UShort_t SubtractPedestal(UShort_t det, 
				    Char_t   rng, 
				    UShort_t sec, 
				    UShort_t str, 
				    Short_t  adc) const;
  /** 
   * Converts number of ADC counts to energy deposited.   This is
   * done by 
   * @f[
   * E_i = A_i g_i
   * @f]
   * where @f$ A_i@f$ is the pedestal subtracted ADC counts, and @f$
   * g_i@f$ is the gain for the @f$ i^{\mbox{th}}@f$ strip. 
   * 
   * @param det	  Detector number  
   * @param rng   Ring identifier 
   * @param sec   Sector number
   * @param str   Strip number 
   * @param eta   Psuedo-rapidity of digit.
   * @param count Pedestal subtracted ADC counts
   *
   * @return Energy deposited @f$ E_i@f$ 
   */
  virtual Float_t  Adc2Energy(UShort_t det, 
			      Char_t   rng, 
			      UShort_t sec, 
			      UShort_t str, 
			      UShort_t count) const;
  /** 
   * Converts number of ADC counts to energy deposited.   This is
   * done by 
   * @f[
   * E_i = A_i g_i
   * @f]
   * where @f$ A_i@f$ is the pedestal subtracted ADC counts, and @f$
   * g_i@f$ is the gain for the @f$ i^{\mbox{th}}@f$ strip. 
   * 
   * @param det	  Detector number  
   * @param rng   Ring identifier 
   * @param sec   Sector number
   * @param str   Strip number 
   * @param eta   Psuedo-rapidity of digit.
   * @param count Pedestal subtracted ADC counts
   *
   * @return Energy deposited @f$ E_i@f$ 
   */
  virtual Float_t  Adc2Energy(UShort_t det, 
			      Char_t   rng, 
			      UShort_t sec, 
			      UShort_t str, 
			      Float_t  eta, 
			      UShort_t count) const;
  /** 
   * Converts an energy signal to number of particles. In this
   * implementation, it's done by 
   * @f[
   * M_i = E_i / E_{\mbox{MIP}}
   * @f]
   * where @f$ E_i@f$ is the energy deposited, and 
   * @f$ E_{\mbox{MIP}}@f$ is the average energy deposited by a
   * minimum ionizing particle
   * 
   * @param det	  Detector number  
   * @param rng   Ring identifier 
   * @param sec   Sector number
   * @param str   Strip number 
   * @param eta   On return, psuedo-rapidity @f$ \eta@f$
   * @param phi   On return, azimuthal angle @f$ \varphi@f$ 
   * @param edep Energy deposited @f$ E_i@f$
   *
   * @return Psuedo-inclusive multiplicity @f$ M@f$ 
   */
  virtual Float_t  Energy2Multiplicity(UShort_t det, 
				       Char_t   rng, 
				       UShort_t sec, 
				       UShort_t str, 
				       Float_t  edep) const;
  /** 
   * Calculate the physical coordinates psuedo-rapidity @f$ \eta@f$,
   * azimuthal angle @f$ \varphi@f$ of the strip corresponding to
   * the digit @a digit.   This is done by using the information
   * obtained, and previously cached by AliFMDGeometry, from the
   * TGeoManager. 
   * 
   * @param det	  Detector number  
   * @param rng   Ring identifier 
   * @param sec   Sector number
   * @param str   Strip number 
   * @param eta   On return, psuedo-rapidity @f$ \eta@f$
   * @param phi   On return, azimuthal angle @f$ \varphi@f$ 
   */
  virtual void     PhysicalCoordinates(UShort_t det, 
				       Char_t   rng, 
				       UShort_t sec, 
				       UShort_t str, 
				       Float_t& eta, 
				       Float_t& phi) const;
  /** 
   * Mark dead channels as invalid, and those that are marked as invalid 
   * but are not dead, get the zero signal. 
   * 
   * @param esd ESD object to modify. 
   */
  void MarkDeadChannels(AliESDFMD* esd) const;

  /** 
   * Utility member function to get the reconstruction parameters for 
   * this event
   * 
   * @return Pointer to AliFMDRecoParam object or null if not
   * available. 
   */
  const AliFMDRecoParam* GetParameters() const;
  /** 
   * Get the numeric identifier of this detector
   * 
   * @return Should be 12
   */  
  Int_t GetIdentifier() const;
  enum Vertex_t {
    kNoVertex,   // Got no vertex
    kGenVertex,  // Got generator vertex 
    kESDVertex   // Got ESD vertex 
  };
  mutable TClonesArray* fMult;          // Cache of RecPoints
  mutable Int_t         fNMult;         // Number of entries in fMult 
  mutable TTree*        fTreeR;         // Output tree 
  mutable Float_t       fCurrentVertex; // Z-coordinate of primary vertex
  mutable AliESDFMD*    fESDObj;        // ESD output object
  mutable Float_t       fNoiseFactor;   // Factor of noise to check
  mutable Bool_t        fAngleCorrect;  // Whether to angle correct
  mutable Vertex_t      fVertexType;    // What kind of vertex we got
  AliESDEvent*          fESD;           // ESD object(?)
  Bool_t                fDiagnostics;   // Wheter to do diagnostics
  TH1*                  fDiagStep1;	// Diagnostics histogram
  TH1*                  fDiagStep2;	// Diagnostics histogram
  TH1*                  fDiagStep3;	// Diagnostics histogram
  TH1*                  fDiagStep4;	// Diagnostics histogram
  TH1*                  fDiagAll;	// Diagnostics histogram
  mutable Bool_t        fZS[3];         // Zero-suppredded?
  mutable UShort_t      fZSFactor[3];   // Noise factor for Zero-suppression
  mutable AliFMDBoolMap fBad;           // Strip marked bad
  Bool_t                fZombie;        // Are we a zombie?
private:
   
  ClassDef(AliFMDReconstructor, 3)  // class for the FMD reconstruction
}; 
#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
// EOF
//
