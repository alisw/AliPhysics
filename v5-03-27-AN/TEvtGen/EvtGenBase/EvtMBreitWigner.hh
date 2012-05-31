#ifndef __EVTMBREITWIGNER_HH__
#define __EVTMBREITWIGNER_HH__

#include "EvtGenBase/EvtMRes.hh"

class EvtMBreitWigner : public EvtMLineShape {

    public:
        EvtMBreitWigner( const EvtId& id, const vector<string>& args );
        ~EvtMBreitWigner( ) {};

        EvtComplex shape( const vector<EvtVector4R>& product ) const;
        
        EvtMLineShape * duplicate() const;

    private:
        EvtId _id;

        double _width, _resmass;
};

#endif
