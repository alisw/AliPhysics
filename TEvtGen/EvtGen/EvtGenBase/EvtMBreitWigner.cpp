#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtMBreitWigner.hh"
#include "EvtGenBase/EvtReport.hh"
#include <stdlib.h>

using std::endl;

EvtMBreitWigner::EvtMBreitWigner( const EvtId& id, const vector<string>& args )
{
    if( args.size() != 0 ) {
        report(Severity::Error, "EvtGen")<<"Unknown input arguments passed in to lineshape."<<endl;
        ::abort();
    }

    _id = id;
    _width = EvtPDL::getWidth( id );
    _resmass = EvtPDL::getMeanMass( id );
}

EvtComplex EvtMBreitWigner::shape( const vector<EvtVector4R>& product ) const 
{
    static EvtComplex I(0.0, 1.0);
    double mass = _node->get4vector(product).mass();
 
    return sqrt(_width/( EvtConst::twoPi )) * 1/( mass - _resmass - I * _width/2 );
}


EvtMLineShape * EvtMBreitWigner::duplicate() const
{
  vector<string> args;
  EvtMLineShape* tmp=new EvtMBreitWigner( _id, args );
  return tmp;
}

