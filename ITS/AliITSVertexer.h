#ifndef ALIITSVERTEXER_H
#define ALIITSVERTEXER_H

#include<AliVertexer.h>
#include "AliITSDetTypeRec.h"

///////////////////////////////////////////////////////////////////
//                                                               //
// Base class for primary vertex reconstruction  for ITS         //
//                                                               //
///////////////////////////////////////////////////////////////////

/* $Id$ */

class TString;

class AliITSVertexer : public AliVertexer {

 public:
    // default constructor
    AliITSVertexer();   
    virtual ~AliITSVertexer();
    virtual AliESDVertex *FindVertexForCurrentEvent(TTree *itsClusterTree)=0;
    virtual void PrintStatus() const = 0;

    virtual void ResetVertex();
    void FindMultiplicity(TTree *itsClusterTree);
    void SetFirstEvent(Int_t ev){fFirstEvent = ev;}
    void SetLastEvent(Int_t ev){fLastEvent = ev;}
    static Float_t GetPipeRadius() {return fgkPipeRadius;}
    void SetLaddersOnLayer2(Int_t ladwid=4);
    void SetComputeMultiplicity(Bool_t opt=kTRUE){fComputeMultiplicity=opt;}
    virtual void SetUseModule(Int_t imod, Bool_t optUse){
      if(imod>=0 && imod<kNSPDMod) fUseModule[imod]=optUse;
    }    
    virtual Bool_t IsModuleUsed(Int_t imod) const {
      if(imod>=0 && imod<kNSPDMod) return fUseModule[imod];
      else return 0;
    }

    virtual AliESDVertex* GetAllVertices(Int_t &novertices) const {novertices = fNoVertices; return fVertArray; }

    AliITSDetTypeRec *GetDetTypeRec() const {return fDetTypeRec;}
    virtual void SetDetTypeRec(AliITSDetTypeRec *ptr){fDetTypeRec = ptr;}
    enum{kNSPDMod=240};

    // Methods containing run-loaders, should be moved to some other class
    void Init(TString filename);
    void WriteCurrentVertex();
    void FindVertices();

    Int_t GetMinTrackletsForPilup() const {return fMinTrackletsForPilup;}
    void SetMinTrackletsForPilup(Int_t mint=2) {fMinTrackletsForPilup=mint;} 
    Bool_t IsPileup() const{return fIsPileup;}
    Float_t GetZPileupV() const {
      if(IsPileup()) return fZpuv;
      else return -9999999.;
    }
    Int_t GetNTrackletsPileupV() const {
      if(IsPileup()) return fNTrpuv;
      else return -2;
    }
 protected:
    static const Float_t fgkPipeRadius;  // beam pipe radius (cm)
    UShort_t *fLadders; // array with layer1-layer2 ladders correspondances  
    Int_t fLadOnLay2;   // (2*fLadOnLay2+1)=number of layer2 ladders 
                      // associated to a layer1 ladder
    Bool_t  fComputeMultiplicity;      // flag to switch on/off tracklet calculation
    Bool_t  fUseModule[kNSPDMod]; // flag for enabling/disabling SPD modules
    AliITSDetTypeRec *fDetTypeRec;  //! pointer to DetTypeRec
    Int_t fMinTrackletsForPilup;  // min. n. of tracklets for pilup definition
    Bool_t fIsPileup;             // flag for pileup
    Int_t fNTrpuv;             // tracklets in pile-up vertex
    Float_t fZpuv;             // Z of second pile-up vertex
    Int_t fNoVertices;         //! number of vertices found 
    AliESDVertex* fVertArray;    //! vertices (main+pileupped)

 private:
    // copy constructor (NO copy allowed: the constructor is protected
    // to avoid misuse)
    AliITSVertexer(const AliITSVertexer& vtxr);
    // assignment operator (NO assignment allowed)
    AliITSVertexer& operator=(const AliITSVertexer& /* vtxr */);

    Int_t fFirstEvent;          // First event to be processed by FindVertices
    Int_t fLastEvent;           // Last event to be processed by FindVertices

  ClassDef(AliITSVertexer,11);
};

#endif
