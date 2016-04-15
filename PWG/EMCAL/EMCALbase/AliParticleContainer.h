#ifndef ALIPARTICLECONTAINER_H
#define ALIPARTICLECONTAINER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <iterator>
#include <TArrayI.h>
#include <AliVParticle.h>

class AliVEvent;
class AliTLorentzVector;

#include "AliEmcalContainer.h"

/**
 * @class AliParticleContainer
 * @brief Container for particles within the EMCAL framework
 * @ingroup EMCALCOREFW
 * @author Martha Verweij
 * @author Salvatore Aiola
 *
 * Container with name, TClonesArray and cuts for particles
 */
class AliParticleContainer : public AliEmcalContainer {
 public:

  /**
   * @class accept_iterator
   * @brief stl iterator over accepted clusters
   * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
   * @date March 16, 2016
   *
   * Bi-directional iterator over all accepted particles in the array. The
   * iterator is used in the standard way a c++ iterator is used. The
   * container is responsible creating it. Afterwards one uses the prefix
   * or postfix increment or decrement operator, depending on the direction
   * of the iteration.
   *
   * ~~~{.cxx}
   * AliParticleContainer cont;
   * for(AliParticleContainer::iterator it = cont.accept_begin(); it != cont.accept_end(); ++cont){
   *   std::cout << "Particle pt: " << (*it)->Pt() std::endl;      // do something with the particle inside
   * }
   * ~~~
   *
   * As the postfix operator needs to perform a copy of the iterator,
   * the prefix operator will have a better performance and should
   * be preferred.
   *
   * Note: This class is using a list of accepted clusters in the
   * backend. This list is copied in the copy constructor and
   * the assignment operator. Can be done in a more performant way
   * using c++11 shared_ptr.
   *
   */
  class accept_iterator : public std::iterator<std::bidirectional_iterator_tag,
                                               AliVParticle,std::ptrdiff_t,
                                               AliVParticle **, AliVParticle *&>
  {
  public:
    accept_iterator(const AliParticleContainer *cont, int startpos, bool forward  = true);
    accept_iterator(const accept_iterator &other);
    virtual ~accept_iterator() {}
    accept_iterator &operator=(const accept_iterator &other);
    bool operator!=(const accept_iterator &other) const;

    accept_iterator &operator++();
    accept_iterator operator++(int);
    accept_iterator &operator--();
    accept_iterator operator--(int);

    AliVParticle *operator*() const;

  private:
    accept_iterator();

    const AliParticleContainer      *fkContainer;       ///< Container iterated over
    TArrayI                          fAcceptIndices;    ///< Indices of accepted clusters
    int                              fCurrentPos;       ///< Current position inside the container
    bool                             fForward;          ///< Direction, expressed in forward direction
  };

  /**
   * @class all_iterator
   * @brief stl iterator over all clusters
   * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
   * @date March 16, 2016
   *
   * Bi-directional iterator over all particles in the array. The
   * iterator is used in the standard way a c++ iterator is used.
   * The container is responsible creating it. Afterwards one uses
   * the prefix or postfix increment or decrement operator, depending
   * on the direction of the iteration.
   *
   * ~~~{.cxx}
   * AliParticleContainer cont;
   * for(AliParticleContainer::iterator it = cont.begin(); it != cont.end(); ++cont){
   *   std::cout << "Particle pt: " << (*it)->Pt() std::endl;      // do something with the particle inside
   * }
   * ~~~
   *
   * If you are using c++11 this reduces to
   *
   * ~~~{.cxx}
   * for(auto it : cont){
   *   std::cout << "Particle pt: " << it->Pt() std::endl;      // do something with the particle inside
   * }
   * ~~~
   *
   * using c++11 range-based iterations.
   *
   * As the postfix operator needs to perform a copy of the iterator,
   * the prefix operator will have a better performance and should
   * be preferred.
   */
  class all_iterator : public std::iterator<std::bidirectional_iterator_tag,
                                            AliVParticle,std::ptrdiff_t,
                                            AliVParticle **, AliVParticle *&>
  {
  public:
    all_iterator(const AliParticleContainer *cont, int startpos, bool forward  = true);
    all_iterator(const all_iterator &other);
    virtual ~all_iterator() {}
    all_iterator &operator=(const all_iterator &other);
    bool operator!=(const all_iterator &other) const;

    all_iterator &operator++();
    all_iterator operator++(int);
    all_iterator &operator--();
    all_iterator operator--(int);

    AliVParticle *operator*() const;

  private:
    all_iterator();

    const AliParticleContainer      *fkContainer;       ///< Container iterated over
    int                              fCurrentPos;       ///< Current position inside the container
    bool                             fForward;          ///< Direction, expressed in forward direction
  };

  AliParticleContainer();
  AliParticleContainer(const char *name);
  virtual ~AliParticleContainer(){;}

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
  void                        SetCharge(Short_t c)                              { fCharge = c         ; }
  void                        SelectHIJING(Bool_t s)                            { if (s) fGeneratorIndex = 0; else fGeneratorIndex = -1; }
  void                        SetGeneratorIndex(Short_t i)                      { fGeneratorIndex = i  ; }

  const char*                 GetTitle() const;

  accept_iterator             accept_begin() const;
  accept_iterator             accept_end() const;
  accept_iterator             accept_rbegin() const;
  accept_iterator             accept_rend() const;

  all_iterator                begin() const;
  all_iterator                end() const;
  all_iterator                rbegin() const;
  all_iterator                rend() const;

 protected:

  Double_t                    fMinDistanceTPCSectorEdge;      ///< require minimum distance to edge of TPC sector edge
  Short_t                     fCharge;                        ///< select particles with charge=fCharge
  Short_t                     fGeneratorIndex;                ///< select MC particles with generator index (default = -1 = switch off selection)

 private:
  AliParticleContainer(const AliParticleContainer& obj); // copy constructor
  AliParticleContainer& operator=(const AliParticleContainer& other); // assignment

  /// \cond CLASSIMP
  ClassDef(AliParticleContainer,9);
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

