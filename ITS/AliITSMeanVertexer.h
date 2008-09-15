#ifndef ALIITSMEANVERTEXER_H
#define ALIITSMEANVERTEXER_H

///////////////////////////////////////////////////////////////////////
//                                                                   //
// Class to compute vertex position using SPD local reconstruction   //
// An average vertex position using all the events                   //
// is built and saved                                                //
///////////////////////////////////////////////////////////////////////

/* $Id$ */

#include <TObject.h>

class TH1F;
class TH2F;
class AliRawReader;
class AliMultiplicity;
class AliESDVertex;
class AliITSDetTypeRec;

class AliITSMeanVertexer : public TObject {

 public:
  // default constructor
    AliITSMeanVertexer();   
    virtual ~AliITSMeanVertexer();
 
    Bool_t Init();
 
    void   SetFilterOnContributors(Int_t nc=1){fFilterOnContributors = nc;}
    void   SetFilterOnTracklets(Int_t nc=1){fFilterOnTracklets = nc;}
    Bool_t Reconstruct(AliRawReader *rawReader, Bool_t mode = kTRUE);
    void   WriteVertices(const char *filename);

    const TH2F*GetVertexXY() const { return fVertexXY; }
    const TH1F*GetVertexZ()  const { return fVertexZ;  }
 
 private:
    // copy constructor (NO copy allowed)
    AliITSMeanVertexer(const AliITSMeanVertexer& vtxr);
    // assignment operator (NO assignment allowed)
    AliITSMeanVertexer& operator=(const AliITSMeanVertexer& /* vtxr */);

    Bool_t Filter(AliESDVertex *vert,AliMultiplicity *mult);
    void   AddToMean(AliESDVertex *vert);
    Bool_t ComputeMean();

    AliITSDetTypeRec *fDetTypeRec; //! ITS reco class
    TH2F    *fVertexXY;        //! histogram with transverse vertex distribution (vertex diamond)
    TH1F    *fVertexZ;         //! histogram with longitudinal vertex distribution

    Double_t fWeighPos[3];     //! weighted average position
    Double_t fWeighSig[3];     //! errors on weighted average position
    Double_t fAverPos[3];      //! average position
    Double_t fAverPosSq[3][3];  //! average square position for covariance
    Int_t fNoEventsContr;      //! number of events used for mean vertex
    Float_t fTotTracklets;       //! total number of tracklets used (integrated)
    Float_t fAverTracklets;    //! average number of tracklets per event
    Float_t fSigmaOnAverTracks; //! RMS of fAverTracklets
    Int_t fFilterOnContributors; //! Numb. of contrib must be > fFilter...
    Int_t fFilterOnTracklets; //! Numb. of tracklets must be > fFilterOnTr...

    ClassDef(AliITSMeanVertexer,0)
};

#endif
