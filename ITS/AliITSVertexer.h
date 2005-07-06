#ifndef ALIITSVERTEXER_H
#define ALIITSVERTEXER_H

#include<AliVertexer.h>

///////////////////////////////////////////////////////////////////
//                                                               //
// Base class for primary vertex reconstruction  for ITS         //
//                                                               //
///////////////////////////////////////////////////////////////////

class TString;
class TClonesArray;


class AliITSVertexer : public AliVertexer {

 public:
    // default constructor
    AliITSVertexer();   
    // standard constructor     
    AliITSVertexer(TString filename); 
    virtual ~AliITSVertexer(){;}
    virtual void SetUseV2Clusters(Bool_t v2c){fUseV2Clusters = v2c;}
    virtual void WriteCurrentVertex();
    virtual void Clusters2RecPoints(const TClonesArray *clusters, Int_t idx, TClonesArray *points);
 

 
 protected:
    // copy constructor (NO copy allowed: the constructor is protected
    // to avoid misuse)
    AliITSVertexer(const AliITSVertexer& vtxr);
    // assignment operator (NO assignment allowed)
    AliITSVertexer& operator=(const AliITSVertexer& /* vtxr */);
    Bool_t fUseV2Clusters;   // true if V2 clusters are used

  ClassDef(AliITSVertexer,2);
};

#endif
