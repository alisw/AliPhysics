#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtResonance.hh"
#include "EvtGenBase/EvtResonance2.hh"
#include "EvtGenModels/EvtMultibody.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtdFunction.hh"
#include "EvtGenBase/EvtKine.hh"
#include "EvtGenBase/EvtParticle.hh"

EvtMultibody::~EvtMultibody() 
{
    if( _decayTree != NULL ) delete _decayTree;
    _decayTree=NULL;
    if( _ilist != NULL ) delete [] _ilist;
    _ilist=NULL;
}

std::string EvtMultibody::getName()
{
    return "D_MULTIBODY";
}

EvtDecayBase* EvtMultibody::clone()
{
    return new EvtMultibody;
}

void EvtMultibody::init()
{
    int N = getNArg();
 
    _decayTree = new EvtMTree( getDaugs(), getNDaug() );
    _ilist = new int[getNDaug()+1];

    for(int i=0; i<N-1; ++i) {
        if(getArgStr( i )=="RESONANCE") {
            _decayTree->addtree( getArgStr( ++i ) );
        } else {
            report(Severity::Error,"EvtGen")
                << "Syntax error at " << getArgStr( i ) << std::endl;
            ::abort();
        }
    }
}

// Set the maximum probability amplitude - if function is left blank then the
// program will search for it.  This however is not deterministic and therefore
// in the release cannot be in place.
void EvtMultibody::initProbMax()
{
    // setProbMax(1.0);
}

void EvtMultibody::decay( EvtParticle *p )
{
    // Initialize the phase space before doing anything else!
    p->initializePhaseSpace(getNDaug(),getDaugs());

    EvtSpinAmp amp = _decayTree->amplitude( p );
    
    vector<int> index = amp.iterallowedinit();
    vector<unsigned int> spins = amp.dims();

    do {
        for( size_t i=0; i<index.size(); ++i ) {
            _ilist[i]=index[i]+spins[i];
        }
        
        vertex( _ilist, amp( index ) );
    } while( amp.iterateallowed( index ) );

}
