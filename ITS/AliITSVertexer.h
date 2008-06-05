#ifndef ALIITSVERTEXER_H
#define ALIITSVERTEXER_H

#include<AliVertexer.h>

///////////////////////////////////////////////////////////////////
//                                                               //
// Base class for primary vertex reconstruction  for ITS         //
//                                                               //
///////////////////////////////////////////////////////////////////

class TString;

class AliITSVertexer : public AliVertexer {

 public:
    // default constructor
    AliITSVertexer();   
    virtual ~AliITSVertexer();
    virtual AliESDVertex *FindVertexForCurrentEvent(TTree *itsClusterTree)=0;
    virtual void PrintStatus() const = 0;

    void FindMultiplicity(TTree *itsClusterTree);
    void SetFirstEvent(Int_t ev){fFirstEvent = ev;}
    void SetLastEvent(Int_t ev){fLastEvent = ev;}
    const Float_t GetPipeRadius()const {return fgkPipeRadius;}
    void SetLaddersOnLayer2(Int_t ladwid=4);

    // Methods containing run-loaders, should be moved to some other class
    void Init(TString filename);
    void WriteCurrentVertex();
    void FindVertices();

 protected:
    static const Float_t fgkPipeRadius;  // beam pipe radius (cm)
    UShort_t *fLadders; // array with layer1-layer2 ladders correspondances  
    Int_t fLadOnLay2;   // (2*fLadOnLay2+1)=number of layer2 ladders 
                      // associated to a layer1 ladder
 
 private:
    // copy constructor (NO copy allowed: the constructor is protected
    // to avoid misuse)
    AliITSVertexer(const AliITSVertexer& vtxr);
    // assignment operator (NO assignment allowed)
    AliITSVertexer& operator=(const AliITSVertexer& /* vtxr */);

    Int_t fFirstEvent;          // First event to be processed by FindVertices
    Int_t fLastEvent;           // Last event to be processed by FindVertices 

  ClassDef(AliITSVertexer,5);
};

#endif
