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

//____________________________________________________________________
class TTree;
class TClonesArray;
class AliFMDDigit;
class AliRawReader;
class AliRunLoader;
class AliESD;
class AliESDFMD;

/** @defgroup FMD_rec Reconstruction */
//____________________________________________________________________
/** @brief This is a class that reconstructs AliFMDRecPoint objects from of
    Digits.  
    This class reads either digits from a TClonesArray or raw data
    from a DDL file (or similar), and applies calibrations to get
    psuedo-inclusive multiplicities per strip.

    @ingroup FMD_rec
 */
class AliFMDReconstructor: public AliReconstructor 
{
public:
  /** CTOR */
  AliFMDReconstructor();
  /** DTOR */
  virtual ~AliFMDReconstructor();

  /** Initialize the reconstructor.  Here, we initialize the geometry
      manager, and finds the local to global transformations from the
      geometry.   The calibration parameter manager is also
      initialized (meaning that the calibration parameters is read
      from CDB).   Next, we try to get some information about the run
      from the run loader passed. 
      @param runLoader Run loader to use to load and store data. 
  */
  virtual void   Init(AliRunLoader* runLoader);
  /** Flag that we can convert raw data into digits. 
      @return always @c true */
  virtual Bool_t HasDigitConversion() const { return kTRUE; }
  /** Convert raw data read from the AliRawReader @a reader into
      digits.  This is done using AliFMDRawReader and
      AliFMDAltroReader.  The digits are put in the passed TTree @a
      digitsTree. 
      @param reader     Raw reader. 
      @param digitsTree Tree to store read digits in. */
  virtual void   ConvertDigits(AliRawReader* reader, TTree* digitsTree) const;
  /** Flag that we can do one-event reconstruction. 
      @return always @c true  */
  virtual Bool_t HasLocalReconstruction() const { return kTRUE; }
  /** Reconstruct one event from the digits passed in @a digitsTree.
      The member function creates AliFMDRecPoint objects and stores
      them on the output tree @a clusterTree.  An FMD ESD object is
      created in parallel. 
      @todo Make sure we get a vertex. 
      @param digitsTree  Tree holding the digits of this event
      @param clusterTree Tree to store AliFMDRecPoint objects in. */
  virtual void   Reconstruct(TTree* digitsTree, TTree* clusterTree) const;
  /** Put in the ESD data, the FMD ESD data.  The object created by
      the Reconstruct member function is copied to the ESD object. 
      @param digitsTree   Tree of digits for this event - not used
      @param clusterTree  Tree of reconstructed points for this event
      - not used. 
      @param esd ESD object to store data in. 
  */
  virtual void   FillESD(TTree* digitsTree, TTree* clusterTree, 
			 AliESD* esd) const;
  /** Not used */
  virtual void   SetESD(AliESD* esd) { fESD = esd; }
  /** Set the noise factor 
      @param f Factor to use */
  virtual void SetNoiseFactor(Float_t f=3) { fNoiseFactor = f; }
  /** Set whether we should do angle correction or nor 
      @param use If true, do angle correction */
  virtual void SetAngleCorrect(Bool_t use=kTRUE) { fAngleCorrect = use; }
protected:
  /** Copy CTOR 
      @param other Object to copy from. */
  AliFMDReconstructor(const AliFMDReconstructor& other);
  /** Assignment operator 
      @param other Object to assign from
      @return reference to this object */
  AliFMDReconstructor& operator=(const AliFMDReconstructor& other);
  /** Process AliFMDDigit objects in @a digits.  For each digit, find
      the psuedo-rapidity @f$ \eta@f$, azimuthal angle @f$ \varphi@f$,
      energy deposited @f$ E@f$, and psuedo-inclusive multiplicity @f$
      M@f$.
      @param digits Array of digits. */
  virtual void     ProcessDigits(TClonesArray* digits) const;
  /** Substract pedestals from raw ADC in @a digit
      @param digit Digit data
      @return Pedestal subtracted ADC count. */
  virtual UShort_t SubtractPedestal(AliFMDDigit* digit) const;
  /** Converts number of ADC counts to energy deposited.   This is
      done by 
      @f[
      E_i = A_i g_i
      @f]
      where @f$ A_i@f$ is the pedestal subtracted ADC counts, and @f$
      g_i@f$ is the gain for the @f$ i^{\mbox{th}}@f$ strip. 
      @param digit Raw data
      @param eta   Psuedo-rapidity of digit.
      @param count Pedestal subtracted ADC counts
      @return Energy deposited @f$ E_i@f$ */
  virtual Float_t  Adc2Energy(AliFMDDigit* digit, Float_t eta, 
			      UShort_t count) const;
  /** Converts an energy signal to number of particles. In this
      implementation, it's done by 
      @f[
      M_i = E_i / E_{\mbox{MIP}}
      @f]
      where @f$ E_i@f$ is the energy deposited, and @f$
      E_{\mbox{MIP}}@f$ is the average energy deposited by a minimum
      ionizing particle
      @param digit Raw data
      @param edep Energy deposited @f$ E_i@f$
      @return Psuedo-inclusive multiplicity @f$ M@f$ */
  virtual Float_t  Energy2Multiplicity(AliFMDDigit* digit, Float_t edep) const;
  /** Calculate the physical coordinates psuedo-rapidity @f$ \eta@f$,
      azimuthal angle @f$ \varphi@f$ of the strip corresponding to
      the digit @a digit.   This is done by using the information
      obtained, and previously cached by AliFMDGeometry, from the
      TGeoManager. 
      @param digit Digit.
      @param eta   On return, psuedo-rapidity @f$ \eta@f$
      @param phi   On return, azimuthal angle @f$ \varphi@f$ */
  virtual void     PhysicalCoordinates(AliFMDDigit* digit, Float_t& eta, 
				       Float_t& phi) const;
  
  mutable TClonesArray* fMult;          // Cache of RecPoints
  mutable Int_t         fNMult;         // Number of entries in fMult 
  mutable TTree*        fTreeR;         // Output tree 
  mutable Float_t       fCurrentVertex; // Z-coordinate of primary vertex
  mutable AliESDFMD*    fESDObj;        // ESD output object
  mutable Float_t       fNoiseFactor;   // Factor of noise to check
  mutable Bool_t        fAngleCorrect;  // Whether to angle correct
  AliESD*               fESD;           // ESD object(?)
  
private:
  /** Hide base classes unused function */
  void Reconstruct(AliRawReader*, TTree*) const;
  /** Hide base classes unused function */
  void Reconstruct(AliRunLoader*) const;
  /** Hide base classes unused function */
  void Reconstruct(AliRunLoader*, AliRawReader*) const;
  /** Hide base classes unused function */
  void FillESD(AliRawReader*, TTree*, AliESD*) const;
  /** Hide base classes unused function */
  void FillESD(AliRunLoader*, AliESD*) const;
  /** Hide base classes unused function */
  void FillESD(AliRunLoader*, AliRawReader*, AliESD*) const;
  
  ClassDef(AliFMDReconstructor, 0)  // class for the FMD reconstruction
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
