#ifndef ALIVERTEXER_H
#define ALIVERTEXER_H

#include<TObject.h>

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
    // computes the vetex for each event and stores it on file
    virtual void FindVertices()= 0;
    virtual void PrintStatus() const = 0;
    virtual void SetDebug(Int_t debug = 0);
    virtual void SetFirstEvent(Int_t ev){fFirstEvent = ev;}
    virtual void SetLastEvent(Int_t ev){fLastEvent = ev;}
    virtual void SetUseV2Clusters(Bool_t choice) = 0;
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
    Int_t fDebug;               //! debug flag - verbose printing if >0

  ClassDef(AliVertexer,1);
};

#endif
