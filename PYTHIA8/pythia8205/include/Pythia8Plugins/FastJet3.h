// FastJet3.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This header file written by Gavin Salam.

#ifndef Pythia8_FastJet3_H
#define Pythia8_FastJet3_H

//----------------------------------------------------------------------
/// \file FastJet3Pythia8.hh
///
/// Code providing an interface for FastJet 3 to make use of Pythia8
/// particles and momenta. Given a
///
/// \code
///   Pythia8::Particle  py8_particle;
/// \endcode
///
/// you may write
///
/// \code
///   fastjet::PseudoJet fj_particle = py8_particle;
/// \endcode
///
/// A copy of the Pythia8::Particle can then be accessed as
///
/// \code
///   fj_particle.user_info<Pythia8::Particle>()
/// \endcode
///
/// so that one can obtain information about the particle such as
///
/// \code
///   fj_particle.user_info<Pythia8::Particle>().status();
///   fj_particle.user_info<Pythia8::Particle>().charge();
/// \endcode
///
/// etc. Note that because the construction of a PseudoJet from the
/// Pythia8 particle involves taking a copy of the whole particle
/// (which has a number of member variables), there will be a small
/// time penalty at that point.
///
/// This file also defines a number of selectors that act on such
/// PseudoJets, such as
///
/// \code
///   SelectorIsCharged();
///   SelectorId(int id);
/// \endcode
///
/// so that one can for example write
///
/// \code
///   vector<PseudoJet> charged_constituents
///     = SelectorIsCharged()(jet.constituents());
/// \endcode
///
/// The full list of Pythia8-specific selectors is to be found at the
/// end of this file. They can be combined with each other and with
/// FastJet selectors using standard boolean operators.  They are all
/// in the fastjet namespace.
///
/// If you do not need the above facilities, then you may instead
/// construct the PseudoJet from the pythia8 particle's 4-vector
///
/// \code
///   PseudoJet fj_particle = py8_particle.p();
/// \endcode
///
/// NB: this code is entirely given as an include file. If compilation
/// time is critical for your application, you may wish to split it
/// into separate .cc and .hh files.
///
// ----------------------------------------------------------------------
// Copyright 2011 by Matteo Cacciari, Gavin Salam and Gregory
// Soyez. Permission is granted to redistribute this file and modify
// it, as long as this notice is retained and any changes are clearly
// marked. No warranties are provided!
// ----------------------------------------------------------------------

#include "fastjet/config.h"             // will allow a test for FJ3
#include "fastjet/ClusterSequence.hh"   // also gives PseudoJet & JetDefinition
#include "fastjet/Selector.hh"
#include "Pythia8/Event.h"              // this is what we need from Pythia8

// FASTJET_VERSION is only defined from version 3 onwards so we can
// use it to test that we have a sufficiently recent version
#ifndef FASTJET_VERSION
#error "FastJet3 is required in order to obtain the features of this interface"
#endif

FASTJET_BEGIN_NAMESPACE // place the code here inside the FJ namespace

/// \class Py8Particle
///
/// A class derived from a pythia 8 particle and that also derives
/// from PseudoJet::UserInfoBase, so that it can be used as UserInfo
/// inside PseudoJets, but also be cast back to the Pythia8 particle
class Py8Particle: public Pythia8::Particle,
                   public PseudoJet::UserInfoBase {
public:
  Py8Particle(const Pythia8::Particle & particle) : Particle(particle) {}
};

/// specialization of the PseudoJet constructor so that it can take a
/// pythia8 particle (and makes a copy of it as user info);
template<>
inline PseudoJet::PseudoJet(const Pythia8::Particle & particle) {
  reset(particle.px(),particle.py(),particle.pz(), particle.e());
  set_user_info(new Py8Particle(particle));
}

/// specialization of the PseudoJet constructor so that it can take a
/// pythia8 Vec4. There is then no particular user info available.
template<>
inline PseudoJet::PseudoJet(const Pythia8::Vec4 & particle) {
  reset(particle.px(),particle.py(),particle.pz(), particle.e());
}


/// \class SelectorWorkerPy8
///
/// A template class to help with the creation of Selectors for Pythia
/// particle properties. It's not necessary to understand how this
/// works in order to use the selectors. See below for the actual list
/// of selectors.
///
/// (But if you're curious, essentially it stores a pointer to a
/// member function of Pythia8::Particle, and when called to select
/// particles, executes it and checks the return value is equal to
/// that requested in the constructor).
template<class T> class SelectorWorkerPy8 : public SelectorWorker {
public:
  /// the typedef helps with the notation for member function pointers
  typedef  T (Pythia8::Particle::*Py8ParticleFnPtr)() const;

  /// c'tor, which takes the member fn pointer and the return value
  /// that it should be equal to
  SelectorWorkerPy8(Py8ParticleFnPtr member_fn_ptr, T value) :
    _member_fn_ptr(member_fn_ptr), _value(value) {};

  /// the one function from SelectorWorker that must be overloaded to
  /// get functioning selection. It makes sure that the PseudoJet
  /// actually has Pythia8::Particle user info before checking
  /// its value.
  bool pass(const PseudoJet & p) const {
    const Pythia8::Particle * py8_particle
      = dynamic_cast<const Pythia8::Particle *>(p.user_info_ptr());
    if (py8_particle == 0) {
      return false; // no info, so false
    } else {
      return (py8_particle->*_member_fn_ptr)() == _value;
    }
  }
private:
  Py8ParticleFnPtr _member_fn_ptr;
  T _value;
};

/// @name Boolean FJ3/PY8 Selectors
///
/// A series of selectors for boolean properties of PseudoJets with
/// Pythia8::Particle information; PseudoJets without
/// Pythia8::Particle structure never pass these selectors.
///
///\{
inline Selector SelectorIsFinal    () {return
  Selector(new SelectorWorkerPy8<bool>(&Pythia8::Particle::isFinal   , true));}
inline Selector SelectorIsCharged  () {return
  Selector(new SelectorWorkerPy8<bool>(&Pythia8::Particle::isCharged , true));}
inline Selector SelectorIsNeutral  () {return
  Selector(new SelectorWorkerPy8<bool>(&Pythia8::Particle::isNeutral , true));}
inline Selector SelectorIsResonance() {return
  Selector(new SelectorWorkerPy8<bool>(&Pythia8::Particle::isResonance,true));}
inline Selector SelectorIsVisible  () {return
  Selector(new SelectorWorkerPy8<bool>(&Pythia8::Particle::isVisible , true));}
inline Selector SelectorIsLepton   () {return
  Selector(new SelectorWorkerPy8<bool>(&Pythia8::Particle::isLepton  , true));}
inline Selector SelectorIsQuark    () {return
  Selector(new SelectorWorkerPy8<bool>(&Pythia8::Particle::isQuark   , true));}
inline Selector SelectorIsGluon    () {return
  Selector(new SelectorWorkerPy8<bool>(&Pythia8::Particle::isGluon   , true));}
inline Selector SelectorIsDiquark  () {return
  Selector(new SelectorWorkerPy8<bool>(&Pythia8::Particle::isDiquark , true));}
inline Selector SelectorIsParton   () {return
  Selector(new SelectorWorkerPy8<bool>(&Pythia8::Particle::isParton  , true));}
inline Selector SelectorIsHadron   () {return
  Selector(new SelectorWorkerPy8<bool>(&Pythia8::Particle::isHadron  , true));}
///\}

/// @name Integer FJ3/PY8 Selectors
///
/// A series of selectors for integer properties of PseudoJets with
/// Pythia8::Particle information; PseudoJets without
/// Pythia8::Particle structure never pass these selectors.
///
///\{
inline Selector SelectorId       (int i) {return
  Selector(new SelectorWorkerPy8<int>(&Pythia8::Particle::id       , i));}
inline Selector SelectorIdAbs    (int i) {return
  Selector(new SelectorWorkerPy8<int>(&Pythia8::Particle::idAbs    , i));}
inline Selector SelectorStatus   (int i) {return
  Selector(new SelectorWorkerPy8<int>(&Pythia8::Particle::status   , i));}
inline Selector SelectorStatusAbs(int i) {return
  Selector(new SelectorWorkerPy8<int>(&Pythia8::Particle::statusAbs, i));}
///\}


FASTJET_END_NAMESPACE

#endif // Pythia8_FastJet3_H
