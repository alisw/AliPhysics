#ifndef EVTMULTIBODY_HH
#define EVTMULTIBODY_HH

#include "EvtGenBase/EvtMTree.hh"
#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtSpinAmp.hh"

class EvtMultibody:public EvtDecayAmp
{
    public:
        EvtMultibody() { _decayTree = NULL; _ilist = NULL; }
        virtual ~EvtMultibody();

        std::string getName();
        EvtDecayBase* clone();

        void init();
        void initProbMax();

        void decay(EvtParticle *p); 

    private:
        EvtMTree * _decayTree;
        int * _ilist;
};

#endif
