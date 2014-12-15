#ifndef __EVTMNODE_HH__
#define __EVTMNODE_HH__

#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtPDL.hh"

#include "EvtGenBase/EvtSymTable.hh"
#include "EvtGenBase/EvtSpinAmp.hh"

#include <vector>
using std::vector;

#include <string>
using std::string;

class EvtMNode {

    public:

        EvtMNode() {}
        virtual ~EvtMNode() {};

        // calculate the amplitude associated event this->children return a
        // vector of the form A_{\lambda this} and sum over allowed angular
        // momenta of the children
        virtual EvtSpinAmp amplitude( const vector<EvtVector4R>
                &product ) const = 0; 

        // get the 4 vector associated with this node
        EvtVector4R get4vector( const vector<EvtVector4R> &product ) const;

        // get twice the spin of the particle
        int getspin() const { return _twospin; }
        EvtSpinType::spintype getspintype() const { return EvtPDL::getSpinType( _id ); }

        // get the id of this node
        EvtId getid() const { return _id; }

        // return which particles this is a combination of
        const vector<int> & getresonance() const { return _resonance; }

        void setparent( EvtMNode * parent ) { _parent = parent; }
        EvtMNode * getparent() const { return _parent; }

        // get the number of children that this node has
        virtual int getnchild() const = 0;
        
        // return the value of the resonance shape
        virtual EvtComplex line( const vector<EvtVector4R>& product ) const=0;
        
        // return a pointer node
        virtual EvtMNode * duplicate() const=0;
    protected:
       
        // store the EvtId of the particle (just in case we need it to access
        // further informatoin about it)
        EvtId _id;
        
        // store TWICE the spin of this resonance (this is to deal with spin 1/2
        int _twospin;

        // store the particles that form this resonance, this should match up
        // with the child nodes from below, and is calculated internally
        vector<int> _resonance;

        // store the parent node of this one
        EvtMNode * _parent;

};

#endif
