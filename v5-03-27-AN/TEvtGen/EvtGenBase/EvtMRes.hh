#ifndef __EVTMRES_HH__
#define __EVTMRES_HH__

#include "EvtGenBase/EvtMNode.hh"

class EvtMRes;

class EvtMLineShape {

    public:

        virtual ~EvtMLineShape() {};

        void setres( EvtMRes * n ) { _node = n; }
        virtual EvtComplex shape( const vector<EvtVector4R>& product ) const=0;
        
        virtual EvtMLineShape * duplicate() const=0;

    protected:

        EvtMRes * _node;
};

class EvtMRes : public EvtMNode {
    
    public:

        virtual ~EvtMRes();

        int getnchild() const { return _children.size(); }

        virtual EvtComplex line( const vector<EvtVector4R>& product ) const 
        { return _lineshape->shape( product ); }

    protected:

        // store the child nodes
        vector<EvtMNode *> _children;

        // store the parametrization amplitudes in some kind 
        EvtSpinAmp _amp;

        // store the lineshape of the resonance
        EvtMLineShape * _lineshape;

};

#endif
