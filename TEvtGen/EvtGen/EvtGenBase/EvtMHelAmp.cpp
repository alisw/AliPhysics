#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtMHelAmp.hh"
#include "EvtGenBase/EvtKine.hh"
#include "EvtGenBase/EvtReport.hh"
#include <stdlib.h>

using std::endl;

EvtMHelAmp::EvtMHelAmp( const EvtId& id, EvtMLineShape * lineshape, 
        const vector<EvtMNode* >& children, const vector<EvtComplex>& elem ) 
{
    _id = id;
    _twospin = EvtSpinType::getSpin2( EvtPDL::getSpinType( id ) );
    _parent = NULL;
    _lineshape = lineshape;
    
    _elem = elem;

    vector<EvtSpinType::spintype> type;
    for(size_t i=0; i<children.size(); ++i ) {
        _children.push_back( children[i] );
        type.push_back( children[i]->getspintype() );
        const vector<int> &res = children[i]->getresonance();
        for(size_t j=0; j<res.size(); ++j )
            _resonance.push_back( res[j] );
        children[i]->setparent( this );
    }
   
    // XXX New code - bugs could appear here XXX
    _amp = EvtSpinAmp( type );
    vector<int> index = _amp.iterinit();
    size_t i = 0;
    do {
        if( !_amp.allowed(index) )
           _amp( index ) = 0.0;
        else if( abs(index[0] - index[1]) > _twospin )
            _amp( index ) = 0.0;
        else {
            _amp( index ) = elem[i];
            ++i;
        }
    } while( _amp.iterate( index ) );
    if(elem.size() != i) {
        report(Severity::Error,"EvtGen")
            <<"Wrong number of elements input in helicity amplitude."<<endl;
        ::abort();
    }

    if( children.size() > 2 ) {
        report(Severity::Error,"EvtGen")
            <<"Helicity amplitude formalism can only handle two body resonances"
            <<endl;
        ::abort();
    }
}

EvtSpinAmp EvtMHelAmp::amplitude( const vector<EvtVector4R> &
        product ) const
{
    EvtVector4R d = _children[0]->get4vector(product);
    double phi, theta;
    
    if( _parent == NULL ) {

        // This means that we're calculating the first level and we need to just
        // calculate the polar and azymuthal angles daughters in rest frame of
        // this (root) particle (this is automatic).
        phi = atan2( d.get(1), d.get(2) ); 
        theta = acos( d.get(3)/d.d3mag() );

    } else {

        // We have parents therefore calculate things in correct coordinate
        // system
        EvtVector4R p = _parent->get4vector(product);
        EvtVector4R q = get4vector(product);

        // See if we have a grandparent - if no then the z-axis is defined by
        // the z-axis of the root particle
        EvtVector4R g = _parent->getparent()==NULL ?
            EvtVector4R(0.0, 0.0, 0.0, 1.0) :
            _parent->getparent()->get4vector(product);

        theta = acos(EvtDecayAngle(p, q, d));
        phi = EvtDecayAnglePhi( g, p, q, d );

    }

    vector<EvtSpinType::spintype> types( 3 );
    types[0] = getspintype();
    types[1] = _children[0]->getspintype();
    types[2] = _children[1]->getspintype();
    EvtSpinAmp amp( types, EvtComplex(0.0, 0.0) );
    vector<int> index = amp.iterallowedinit();

    do {
        if( abs(index[1]-index[2]) > _twospin ) continue;
        amp(index) +=
            conj(wignerD(_twospin,index[0],index[1]-index[2],phi,theta,0.0)) *
            _amp(index[1],index[2]); 
    } while(amp.iterateallowed(index));

    EvtSpinAmp amp0 = _children[0]->amplitude(product);
    EvtSpinAmp amp1 = _children[1]->amplitude(product);

    amp.extcont( amp0, 1, 0 );
    amp.extcont( amp1, 1, 0 );

    amp *= sqrt( ( _twospin + 1 ) / ( 2 * EvtConst::twoPi ) ) *
        _children[0]->line(product) * _children[1]->line(product);

    return amp;
}

EvtMNode * EvtMHelAmp::duplicate() const
{
    vector<EvtMNode *> children;

    for(size_t i=0; i<_children.size(); ++i ) {
        children.push_back( _children[i]->duplicate() );
    }
    
    EvtMLineShape * lineshape = _lineshape->duplicate();
    EvtMHelAmp * ret = new EvtMHelAmp( _id, lineshape, children, _elem );
    lineshape->setres( ret );

    return ret;
}
