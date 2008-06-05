#ifndef ALIITSMEANVERTEXER_H
#define ALIITSMEANVERTEXER_H

#include <TString.h>

///////////////////////////////////////////////////////////////////////
//                                                                   //
// Class to compute vertex position using SPD local reconstruction   //
// An average vertex position using all the events                   //
// is built and saved                                                //
///////////////////////////////////////////////////////////////////////

/* $Id$ */

class TObject;
class AliRawReader;
class AliRunLoader;
class AliMultiplicity;
class AliESDVertex;

class AliITSMeanVertexer : public TObject {

 public:
    // default constructor
    AliITSMeanVertexer();   
    // standard constructor. filename is the name of the file containing
    // raw data either in ROOT or DATE format according to file extension
    AliITSMeanVertexer(TString &filename);   
    // constructor with explicit assignment of names for geometry and loaders
    AliITSMeanVertexer(TString &filename, TString &loaderfilename, 
                       TString &geometryfilename);
    virtual ~AliITSMeanVertexer();
    void SetLoaderFileName(TString fn="ITSMeanVertexer.root")
                           {fLoaderFileName = fn;}
    void SetGeometryFileName(TString fn="geometry.root")
                           {fGeometryFileName = fn;}
    void SetMeanVertexFileName(TString fn) {fMVFileName = fn;}
    void SetFilterOnContributors(Int_t nc=1){fFilterOnContributors = nc;}
    void SetFilterOnTracklets(Int_t nc=1){fFilterOnTracklets = nc;}
    void SetWriteVertices(Bool_t action){fWriteVertices = action;}
    Bool_t GetWriteVertices() const {return fWriteVertices;}
    void Reconstruct();
    void DoVertices();

 
 protected:
    // copy constructor (NO copy allowed: the constructor is protected
    // to avoid misuse)
    AliITSMeanVertexer(const AliITSMeanVertexer& vtxr);
    // assignment operator (NO assignment allowed)
    AliITSMeanVertexer& operator=(const AliITSMeanVertexer& /* vtxr */);
    void Init(TString &filename);  // initialization invoked by constructors
    Bool_t Filter(AliESDVertex *vert,AliMultiplicity *mult);
    void AddToMean(AliESDVertex *vert);
    Bool_t ComputeMean();

    static const TString fgkMVFileNameDefault;  //! default for fMVFileName
    TString fLoaderFileName;   //! name of the local file containing loaders
    TString fGeometryFileName; //! name of the file containing the geometry
    TString fMVFileName;       //! name of the file containing the mean vertex
    AliRawReader *fRawReader;  //! rawreader object
    AliRunLoader *fRunLoader;  //! run loader
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
    Bool_t fWriteVertices; //! if kTRUE all the vertices are saved to a file

  ClassDef(AliITSMeanVertexer,0);
};

#endif
