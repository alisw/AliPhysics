#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtMNode.hh"

EvtVector4R EvtMNode::get4vector( const vector<EvtVector4R> &product )
    const 
{

    EvtVector4R res(0.0, 0.0, 0.0, 0.0);
    vector<int>::const_iterator iter;
    
    for( iter = _resonance.begin(); iter != _resonance.end(); ++iter )
        res += product[ *iter ];
    
    return res;

}
