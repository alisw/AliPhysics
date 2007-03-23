#ifndef ALIVERTEXER_H
#define ALIVERTEXER_H

#include<TObject.h>
#include<AliMultiplicity.h>

///////////////////////////////////////////////////////////////////
//                                                               //
// Base class for primary vertex reconstruction                  //
//                                                               //
///////////////////////////////////////////////////////////////////

class TFile;
class TString;
class TTRee;
class AliESDVertex;


class AliVertexer : public TObject {

 public:
    // default constructor
    AliVertexer();  
 
    // destructor
    virtual ~AliVertexer(); 
    // computes the vertex for the current event
    virtual AliESDVertex* FindVertexForCurrentEvent(Int_t evnumb)=0; 
    // computes the vertex for each event and stores it on file
    virtual void FindVertices()= 0;
    virtual AliMultiplicity* GetMultiplicity() const {return fMult;}
    virtual void PrintStatus() const = 0;
    virtual void SetFirstEvent(Int_t ev){fFirstEvent = ev;}
    virtual void SetLastEvent(Int_t ev){fLastEvent = ev;}
    virtual void SetVtxStart(Double_t x,Double_t y,Double_t z) 
           { fNominalPos[0]=x; fNominalPos[1]=y; fNominalPos[2]=z; }
    virtual void SetVtxStartSigma(Double_t sx,Double_t sy,Double_t sz) 
           { fNominalCov[0]=sx*sx; fNominalCov[2]=sy*sy; fNominalCov[5]=sz*sz;
             fNominalCov[1]=0.; fNominalCov[3]=0.; fNominalCov[4]=0.; }
    virtual void SetVtxStart(AliESDVertex *vtx);
    virtual void WriteCurrentVertex() = 0;

 
 protected:
    // copy constructor (NO copy allowed: the constructor is protected
    // to avoid misuse)
    AliVertexer(const AliVertexer& vtxr);
    // assignment operator (NO assignment allowed)
    AliVertexer& operator=(const AliVertexer& /* vtxr */);

    AliESDVertex *fCurrentVertex;  //! pointer to the vertex of the current
                                   //  event
    Int_t fFirstEvent;          // First event to be processed by FindVertices
    Int_t fLastEvent;           // Last event to be processed by FindVertices 
    Double_t  fNominalPos[3];   // initial knowledge on vertex position
    Double_t  fNominalCov[6];   // initial knowledge on vertex position
    AliMultiplicity *fMult;     //! Multiplicity object

  ClassDef(AliVertexer,3);
};

#endif
