#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtMRes.hh"

EvtMRes::~EvtMRes()
{
    for(size_t i=0; i<_children.size(); ++i)
        delete _children[i];
}
