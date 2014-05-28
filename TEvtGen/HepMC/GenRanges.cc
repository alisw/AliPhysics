//--------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////
// garren@fnal.gov, March 2010
// 
//////////////////////////////////////////////////////////////////////////
//--------------------------------------------------------------------------

#include <iostream>

#include "HepMC/GenRanges.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"

namespace HepMC {

GenEventVertexRange GenEvent::vertex_range()
{
    return GenEventVertexRange(*this);
}

ConstGenEventVertexRange GenEvent::vertex_range() const
{
    return ConstGenEventVertexRange(*this);
}

GenEventParticleRange GenEvent::particle_range()
{
    return GenEventParticleRange(*this);
}

ConstGenEventParticleRange GenEvent::particle_range() const
{
    return ConstGenEventParticleRange(*this);
}

GenVertexParticleRange GenVertex::particles( IteratorRange range )
{
    return GenVertexParticleRange(*this,range);
}

GenParticleProductionRange GenVertex::particles_in( GenParticle& p, IteratorRange range )
{
    return GenParticleProductionRange(p,range);
}

ConstGenParticleProductionRange GenVertex:: particles_in( GenParticle const & p, IteratorRange range ) const
{
    return ConstGenParticleProductionRange(p,range);
}

GenParticleEndRange GenVertex::particles_out( GenParticle& p, IteratorRange range )
{
    return GenParticleEndRange(p,range);
}

ConstGenParticleEndRange GenVertex::particles_out( GenParticle const & p, IteratorRange range ) const
{
    return ConstGenParticleEndRange(p,range);
}

GenParticleProductionRange GenParticle::particles_in( IteratorRange range )
{
    return GenParticleProductionRange(*this,range);
}


ConstGenParticleProductionRange GenParticle::particles_in( IteratorRange range ) const
{
    return ConstGenParticleProductionRange(*this,range);
}


GenParticleEndRange GenParticle::particles_out( IteratorRange range )
{
    return GenParticleEndRange(*this,range);
}


ConstGenParticleEndRange GenParticle::particles_out( IteratorRange range ) const
{
    return ConstGenParticleEndRange(*this,range);
}



} // HepMC
