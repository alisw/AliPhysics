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
    virtual ~AliITSVertexer();
    virtual void FindMultiplicity(Int_t evnumber);
    virtual void WriteCurrentVertex();
    const Float_t GetPipeRadius()const {return fgkPipeRadius;}
    virtual void SetLaddersOnLayer2(Int_t ladwid=4);

 
 protected:
    // copy constructor (NO copy allowed: the constructor is protected
    // to avoid misuse)
    AliITSVertexer(const AliITSVertexer& vtxr);
    // assignment operator (NO assignment allowed)
    AliITSVertexer& operator=(const AliITSVertexer& /* vtxr */);

    static const Float_t fgkPipeRadius;  // beam pipe radius (cm)
    UShort_t *fLadders; // array with layer1-layer2 ladders correspondances  
    Int_t fLadOnLay2;   // (2*fLadOnLay2+1)=number of layer2 ladders 
                      // associated to a layer1 ladder

  ClassDef(AliITSVertexer,4);
};

#endif
