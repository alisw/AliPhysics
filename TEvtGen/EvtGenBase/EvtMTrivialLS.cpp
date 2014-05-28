#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtMTrivialLS.hh"

EvtComplex EvtMTrivialLS::shape( const vector<EvtVector4R>& /*product*/ ) const 
{ 
    return EvtComplex(1.0, 0.0); 
}

EvtMLineShape* EvtMTrivialLS::duplicate() const 
{   
    EvtId temp1;
    vector<string> temp2;

    EvtMLineShape* tmp=new EvtMTrivialLS(temp1, temp2); 
    return tmp;
}
