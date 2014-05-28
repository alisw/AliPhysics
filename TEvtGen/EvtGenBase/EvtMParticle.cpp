#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtMParticle.hh"
#include "EvtGenBase/EvtSpinType.hh"

EvtMParticle::EvtMParticle( int label, const EvtId& id )
{
    _id = id;
    _twospin = EvtSpinType::getSpin2( EvtPDL::getSpinType( id ) );
    _resonance.push_back( label );
}

EvtSpinAmp EvtMParticle::amplitude( const vector<EvtVector4R> &/*product*/) const
{
    vector<EvtSpinType::spintype> types( 2, getspintype() );
    EvtSpinAmp amp( types, EvtComplex( 0.0, 0.0 ) );

    for( int i=-_twospin; i<=_twospin; i+=2 )
        amp(i, i) = EvtComplex( 1.0, 0.0 );

    return amp;
}

EvtMNode * EvtMParticle::duplicate() const
{
    return new EvtMParticle( _resonance[0], _id );
}
