#ifndef ALIPARTICLECONTAINER_H
#define ALIPARTICLECONTAINER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TArrayI.h>
#include <AliVParticle.h>

class AliVEvent;
class AliTLorentzVector;

#include "AliEmcalContainer.h"
#if !(defined(__CINT__) || defined(__MAKECINT__))
#include "AliEmcalContainerIndexMap.h"
#endif

#if !(defined(__CINT__) || defined(__MAKECINT__))
typedef EMCALIterableContainer::AliEmcalIterableContainerT<AliVParticle, EMCALIterableContainer::operator_star_object<AliVParticle> > AliParticleIterableContainer;
typedef EMCALIterableContainer::AliEmcalIterableContainerT<AliVParticle, EMCALIterableContainer::operator_star_pair<AliVParticle> > AliParticleIterableMomentumContainer;
#endif

/**
 * @class AliParticleContainer
 * @brief Container for particles within the EMCAL framework
 * @ingroup EMCALCOREFW
 * @author Marta Verweij
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 *
 * Container with name, TClonesArray and cuts for particles
 */
class AliParticleContainer : public AliEmcalContainer {
 public:
  enum EChargeCut_t { kNoChargeCut, kCharged, kNeutral, kPositiveCharge, kNegativeCharge };

  AliParticleContainer();
  AliParticleContainer(const char *name);
  virtual ~AliParticleContainer(){;}

  /**
   * Index operator: Providing access to track in the container with the
   * given index.
   * @param[in] index Index of the particle in the array
   * @return Particle at the given index
   */
  virtual TObject *operator[](int index) const { return GetParticle(index); }

  virtual Bool_t              ApplyParticleCuts(const AliVParticle* vp, UInt_t &rejectionReason) const;
  virtual Bool_t              ApplyKinematicCuts(const AliTLorentzVector& mom, UInt_t &rejectionReason) const;
  virtual Bool_t              AcceptObject(Int_t i, UInt_t &rejectionReason) const              { return AcceptParticle(i, rejectionReason);}
  virtual Bool_t              AcceptObject(const TObject* obj, UInt_t &rejectionReason) const   { return AcceptParticle(dynamic_cast<const AliVParticle*>(obj), rejectionReason);}
  virtual Bool_t              AcceptParticle(const AliVParticle* vp, UInt_t &rejectionReason) const        ;
  virtual Bool_t              AcceptParticle(Int_t i, UInt_t &rejectionReason) const                       ;
  Double_t                    GetParticlePtCut()                        const   { return GetMinPt()     ; }
  Double_t                    GetParticleEtaMin()                       const   { return GetMinEta()    ; }
  Double_t                    GetParticleEtaMax()                       const   { return GetMaxEta()    ; }
  Double_t                    GetParticlePhiMin()                       const   { return GetMinPhi()    ; }
  Double_t                    GetParticlePhiMax()                       const   { return GetMaxPhi()    ; }
  void                        SetParticlePtCut(Double_t cut)                    { SetMinPt(cut)         ; }
  void                        SetParticleEtaLimits(Double_t min, Double_t max)  { SetEtaLimits(min, max); }
  void                        SetParticlePhiLimits(Double_t min, Double_t max)  { SetPhiLimits(min, max); }
  virtual AliVParticle       *GetLeadingParticle(const char* opt="")         ;
  virtual AliVParticle       *GetParticle(Int_t i=-1)                   const;
  virtual AliVParticle       *GetAcceptParticle(Int_t i=-1)             const;
  virtual AliVParticle       *GetNextAcceptParticle()                        ;
  virtual AliVParticle       *GetNextParticle()                              ;
  virtual Bool_t              GetMomentumFromParticle(TLorentzVector &mom, const AliVParticle* part, Double_t mass) const;
  virtual Bool_t              GetMomentumFromParticle(TLorentzVector &mom, const AliVParticle* part) const;
  virtual Bool_t              GetMomentum(TLorentzVector &mom, Int_t i) const;
  virtual Bool_t              GetAcceptMomentum(TLorentzVector &mom, Int_t i) const;
  virtual Bool_t              GetNextMomentum(TLorentzVector &mom);
  virtual Bool_t              GetNextAcceptMomentum(TLorentzVector &mom);
  Int_t                       GetNParticles()                           const   {return GetNEntries();}
  Int_t                       GetNAcceptedParticles()                   const;
  void                        SetMinDistanceTPCSectorEdge(Double_t min)         { fMinDistanceTPCSectorEdge = min; }
  void                        SetCharge(EChargeCut_t c)                         { fChargeCut = c       ; }
  void                        SelectHIJING(Bool_t s)                            { if (s) fGeneratorIndex = 0; else fGeneratorIndex = -1; }
  void                        SetGeneratorIndex(Short_t i)                      { fGeneratorIndex = i  ; }
  void                        SetArray(const AliVEvent * event);

  const char*                 GetTitle() const;

#if !(defined(__CINT__) || defined(__MAKECINT__))
  /// Get the EMCal container utils associated with particle containers
  static const AliEmcalContainerIndexMap <TClonesArray, AliVParticle>& GetEmcalContainerIndexMap() { return fgEmcalContainerIndexMap; }

  const AliParticleIterableContainer      all() const;
  const AliParticleIterableContainer      accepted() const;

  const AliParticleIterableMomentumContainer      all_momentum() const;
  const AliParticleIterableMomentumContainer      accepted_momentum() const;
#endif

 protected:

#if !(defined(__CINT__) || defined(__MAKECINT__))
  static AliEmcalContainerIndexMap <TClonesArray, AliVParticle> fgEmcalContainerIndexMap; //!<! Mapping from containers to indices
#endif

  Double_t                    fMinDistanceTPCSectorEdge;      ///< require minimum distance to edge of TPC sector edge
  EChargeCut_t                fChargeCut;                     ///< select particles according to their charge
  Short_t                     fGeneratorIndex;                ///< select MC particles with generator index (default = -1 = switch off selection)

 private:
  AliParticleContainer(const AliParticleContainer& obj); // copy constructor
  AliParticleContainer& operator=(const AliParticleContainer& other); // assignment

  /// \cond CLASSIMP
  ClassDef(AliParticleContainer,11);
  /// \endcond

};

/**
 * Unit test for the iterators. Comparing iterators against for-loop of particles.
 * All particles selected in the for-loop must be found in order to pass the test.
 * @ingroup EMCALCOREFW
 * @param cont Particle container used for the test.
 * @param iteratorType type of the iterator (0 = accept_iterator, 1 = all_iterator)
 * @param verbose Switch on verbosity in case of true
 * @return Result of the unit test (0 - passed, 1 - particles missing, 2 - excess particles)
 */
int TestParticleContainerIterator(const AliParticleContainer *const cont, int iteratorType = 0, bool verbose = false);

#endif

