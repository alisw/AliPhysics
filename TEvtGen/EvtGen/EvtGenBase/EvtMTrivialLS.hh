#ifndef __EVTMTRIVIALLS_H__
#define __EVTMTRIVIALLS_H__

#include "EvtGenBase/EvtMRes.hh"

class EvtMTrivialLS : public EvtMLineShape {

    public:

        EvtMTrivialLS( const EvtId& /*id*/, const vector<string>& /*args*/ ) {};
        ~EvtMTrivialLS() {};

        EvtComplex shape( const vector<EvtVector4R>& product ) const; 
    
        EvtMLineShape* duplicate() const;
};

#endif
