#ifndef HEPMC_GEN_EVENT_ITERATORS_H
#define HEPMC_GEN_EVENT_ITERATORS_H

//--------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////
// garren@fnal.gov, May 2009
// 
//////////////////////////////////////////////////////////////////////////
//--------------------------------------------------------------------------

#include <stdexcept>

#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"

namespace HepMC {

//! GenEventVertexRange acts like a collection of vertices

///
/// \class  GenEventVertexRange
/// HepMC::GenEventVertexRange is used to mimic a collection of
/// vertices for ease of use - especially with utilities such as 
/// the Boost foreach funtion
///
class GenEventVertexRange {

public:

  /// the constructor requires a GenEvent
  GenEventVertexRange( GenEvent & e ) : m_event(e) {}
  /// 
  GenEvent::vertex_iterator begin() { return m_event.vertices_begin(); }
  GenEvent::vertex_iterator end()   { return m_event.vertices_end(); }

private:
  /// Because the class contains a reference, assignments are not allowed.
  /// However, we need the copy constructor for GenEvent::vertex_range().
  GenEventVertexRange& operator=( GenEventVertexRange & );

private:
  GenEvent & m_event;

};

//! ConstGenEventVertexRange acts like a collection of vertices

///
/// \class  ConstGenEventVertexRange
/// HepMC::ConstGenEventVertexRange is used to mimic a collection of
/// vertices for ease of use - especially with utilities such as 
/// the Boost foreach funtion
/// This is the const partner of GenEventVertexRange
///
class ConstGenEventVertexRange {

public:

  /// the constructor requires a const GenEvent
  ConstGenEventVertexRange( GenEvent const & e ) : m_event(e) {}
  /// 
  GenEvent::vertex_const_iterator begin() const { return m_event.vertices_begin(); }
  GenEvent::vertex_const_iterator end()   const { return m_event.vertices_end(); }

private:
  /// Because the class contains a reference, assignments are not allowed.
  /// However, we need the copy constructor for GenEvent::vertex_range().
  ConstGenEventVertexRange& operator=( ConstGenEventVertexRange & );

private:
  GenEvent const & m_event;

};

//! GenEventParticleRange acts like a collection of particles

///
/// \class  GenEventParticleRange
/// HepMC::GenEventParticleRange is used to mimic a collection of
/// particles for ease of use - especially with utilities such as 
/// the Boost foreach funtion
///
class GenEventParticleRange {

public:

  /// the constructor requires a GenEvent
  GenEventParticleRange( GenEvent & e ) : m_event(e) {}
  /// 
  GenEvent::particle_iterator begin() { return m_event.particles_begin(); }
  GenEvent::particle_iterator end()   { return m_event.particles_end(); }

private:
  /// Because the class contains a reference, assignments are not allowed.
  /// However, we need the copy constructor for GenEvent::particle_range().
  GenEventParticleRange& operator=( GenEventParticleRange & );

private:
  GenEvent & m_event;

};

//! ConstGenEventParticleRange acts like a collection of particles

///
/// \class  ConstGenEventParticleRange
/// HepMC::ConstGenEventParticleRange is used to mimic a collection of
/// particles for ease of use - especially with utilities such as 
/// the Boost foreach funtion
/// This is the const partner of GenEventParticleRange
///
class ConstGenEventParticleRange {

public:

  /// the constructor requires a const GenEvent
  ConstGenEventParticleRange( GenEvent const & e ) : m_event(e) {}
  /// 
  GenEvent::particle_const_iterator begin() const { return m_event.particles_begin(); }
  GenEvent::particle_const_iterator end()   const { return m_event.particles_end(); }

private:
  /// Because the class contains a reference, assignments are not allowed.
  /// However, we need the copy constructor for GenEvent::particle_range().
  ConstGenEventParticleRange& operator=( ConstGenEventParticleRange & );

private:
  GenEvent const & m_event;

};

//! GenVertexParticleRange acts like a collection of particles

///
/// \class  GenVertexParticleRange
/// HepMC::GenVertexParticleRange is used to mimic a collection of
/// particles for ease of use - especially with utilities such as 
/// the Boost foreach funtion
///
class GenVertexParticleRange {

public:

  /// the constructor requires a GenVertex
  GenVertexParticleRange( GenVertex & v, IteratorRange range = relatives ) 
  : m_vertex(v),m_range(range) {}
  /// 
  GenVertex::particle_iterator begin() { return m_vertex.particles_begin(m_range); }
  GenVertex::particle_iterator end()   { return m_vertex.particles_end(m_range); }

private:
  /// Because the class contains a reference, assignments are not allowed.
  /// However, we need the copy constructor for GenVertex::particles().
  GenVertexParticleRange& operator=( GenVertexParticleRange & );

private:
  GenVertex     & m_vertex;
  IteratorRange   m_range;

};

//! GenParticleProductionRange acts like a collection of particles

///
/// \class  GenParticleProductionRange
/// HepMC::GenParticleProductionRange is used to mimic a collection of
/// particles associated with the particle's production vertex for ease of use 
/// Utilities such as the Boost foreach funtion will want to use this class.
///
class GenParticleProductionRange {

public:

  /// the constructor requires a GenParticle
  GenParticleProductionRange( GenParticle const & p, IteratorRange range = relatives ) 
  : m_particle(p),m_range(range) {}
  /// begin iterator throws an error if the particle production_vertex is undefined
  GenVertex::particle_iterator begin();
  /// end iterator throws an error if the particle production_vertex is undefined
  GenVertex::particle_iterator end(); 

private:
  /// Because the class contains a reference, assignments are not allowed.
  /// However, we need the copy constructor for GenVertex::particles_in().
  GenParticleProductionRange& operator=( GenParticleProductionRange & );

private:
  GenParticle const & m_particle;
  IteratorRange       m_range;

};

class ConstGenParticleProductionRange {

public:

  /// the constructor requires a GenParticle
  ConstGenParticleProductionRange( GenParticle const & p, IteratorRange range = relatives ) 
  : m_particle(p),m_range(range) {}
  /// begin iterator throws an error if the particle production_vertex is undefined
  GenVertex::particle_iterator begin();
  /// end iterator throws an error if the particle production_vertex is undefined
  GenVertex::particle_iterator end(); 

private:
  /// Because the class contains a reference, assignments are not allowed.
  /// However, we need the copy constructor for GenVertex::particles_in().
  ConstGenParticleProductionRange& operator=( ConstGenParticleProductionRange & );

private:
  GenParticle const & m_particle;
  IteratorRange       m_range;

};

//! GenParticleEndRange acts like a collection of particles

///
/// \class  GenParticleEndRange
/// HepMC::GenParticleEndRange is used to mimic a collection of
/// particles associated with the particle's end vertex for ease of use 
/// Utilities such as the Boost foreach funtion will want to use this class.
///
class GenParticleEndRange {

public:

  /// the constructor requires a GenParticle
  GenParticleEndRange( GenParticle const & p, IteratorRange range = relatives ) 
  : m_particle(p),m_range(range) {}
  /// begin iterator throws an error if the particle end_vertex is undefined
  GenVertex::particle_iterator begin();
  /// end iterator throws an error if the particle end_vertex is undefined
  GenVertex::particle_iterator end(); 

private:
  /// Because the class contains a reference, assignments are not allowed.
  /// However, we need the copy constructor for GenVertex::particles_out().
  GenParticleEndRange& operator=( GenParticleEndRange & );

private:
  GenParticle const & m_particle;
  IteratorRange       m_range;

};

class ConstGenParticleEndRange {

public:

  /// the constructor requires a GenParticle
  ConstGenParticleEndRange( GenParticle const & p, IteratorRange range = relatives ) 
  : m_particle(p),m_range(range) {}
  /// begin iterator throws an error if the particle end_vertex is undefined
  GenVertex::particle_iterator begin();
  /// end iterator throws an error if the particle end_vertex is undefined
  GenVertex::particle_iterator end(); 

private:
  /// Because the class contains a reference, assignments are not allowed.
  /// However, we need the copy constructor for GenVertex::particles_out().
  ConstGenParticleEndRange& operator=( ConstGenParticleEndRange & );

private:
  GenParticle const & m_particle;
  IteratorRange       m_range;

};


inline GenVertex::particle_iterator GenParticleProductionRange::begin() 
{ 
    if ( ! m_particle.production_vertex() ) 
        throw(std::range_error("GenParticleProductionRange: GenParticle has no production_vertex"));
    return m_particle.production_vertex()->particles_begin(m_range); 
}

inline GenVertex::particle_iterator GenParticleProductionRange::end()
{ 
    if ( ! m_particle.production_vertex() ) 
        throw(std::range_error("GenParticleProductionRange: GenParticle has no production_vertex"));
    return m_particle.production_vertex()->particles_end(m_range); 
}


inline GenVertex::particle_iterator ConstGenParticleProductionRange::begin() 
{ 
    if ( ! m_particle.production_vertex() ) 
        throw(std::range_error("ConstGenParticleProductionRange: GenParticle has no production_vertex"));
    return m_particle.production_vertex()->particles_begin(m_range); 
}

inline GenVertex::particle_iterator ConstGenParticleProductionRange::end()
{ 
    if ( ! m_particle.production_vertex() ) 
        throw(std::range_error("ConstGenParticleProductionRange: GenParticle has no production_vertex"));
    return m_particle.production_vertex()->particles_end(m_range); 
}

inline GenVertex::particle_iterator GenParticleEndRange::begin() 
{ 
    if ( ! m_particle.end_vertex() ) 
        throw(std::range_error("GenParticleEndRange: GenParticle has no end_vertex"));
    return m_particle.end_vertex()->particles_begin(m_range); 
}
inline GenVertex::particle_iterator GenParticleEndRange::end()
{ 
    if ( ! m_particle.end_vertex() ) 
        throw(std::range_error("GenParticleEndRange: GenParticle has no end_vertex"));
    return m_particle.end_vertex()->particles_end(m_range); 
}

inline GenVertex::particle_iterator ConstGenParticleEndRange::begin() 
{ 
    if ( ! m_particle.end_vertex() ) 
        throw(std::range_error("ConstGenParticleEndRange: GenParticle has no end_vertex"));
    return m_particle.end_vertex()->particles_begin(m_range); 
}
inline GenVertex::particle_iterator ConstGenParticleEndRange::end()
{ 
    if ( ! m_particle.end_vertex() ) 
        throw(std::range_error("ConstGenParticleEndRange: GenParticle has no end_vertex"));
    return m_particle.end_vertex()->particles_end(m_range); 
}

} // HepMC

#endif  // HEPMC_GEN_EVENT_ITERATORS_H
