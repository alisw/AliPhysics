#ifndef ALIVERTEXER_H
#define ALIVERTEXER_H

///////////////////////////////////////////////////////////////////
//                                                               //
// Base class for primary vertex reconstruction                  //
//                                                               //
///////////////////////////////////////////////////////////////////

#include<TObject.h>

class TTree;
class AliESDVertex;
class AliMultiplicity;


class AliVertexer : public TObject {

 public:
    // default constructor
    AliVertexer();  
 
    // destructor
    virtual ~AliVertexer();
    // computes the vertex for the current event
    virtual AliESDVertex* FindVertexForCurrentEvent(TTree *clustersTree)=0; 
    // computes the vertex for each event and stores it on file
    AliMultiplicity* GetMultiplicity() const {return fMult;}
    virtual void PrintStatus() const = 0;
    void SetVtxStart(Double_t x,Double_t y,Double_t z) 
      { fNominalPos[0]=x; fNominalPos[1]=y; fNominalPos[2]=z; }
    void SetVtxStartSigma(Double_t sx,Double_t sy,Double_t sz) 
      { fNominalCov[0]=sx*sx; fNominalCov[2]=sy*sy; fNominalCov[5]=sz*sz;
      fNominalCov[1]=0.; fNominalCov[3]=0.; fNominalCov[4]=0.; }
    void SetVtxStart(AliESDVertex *vtx);
    // the following method can be implemented in daughter classes 
    // (e.g. in AliITSVertexer3D). It is intended to tag pile-up events
    // novertices is the total number of vertices (1 means no pileup)
    // The returned pointer points to an array of AliESDVertx opbjects
    // with size=novertices
    virtual  AliESDVertex* GetAllVertices(Int_t &novertices) const {novertices = 0; return NULL;}
    const Double_t* GetNominalPos() const {return fNominalPos;}
    const Double_t* GetNominalCov() const {return fNominalCov;}

 protected:
    AliESDVertex *fCurrentVertex; //! pointer to the current vertex
    AliMultiplicity *fMult;     //! Multiplicity object
 
 private:
    // copy constructor (NO copy allowed: the constructor is protected
    // to avoid misuse)
    AliVertexer(const AliVertexer& vtxr);
    // assignment operator (NO assignment allowed)
    AliVertexer& operator=(const AliVertexer& /* vtxr */);

    Double_t  fNominalPos[3];   // initial knowledge on vertex position
    Double_t  fNominalCov[6];   // initial knowledge on vertex position

  ClassDef(AliVertexer,4);
};

#endif
