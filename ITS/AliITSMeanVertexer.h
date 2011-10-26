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
#include <TBits.h>

class TClonesArray;
class TH1F;
class TH2F;
class AliRawReader;
class AliESDVertex;
class AliITSDetTypeRec;
class AliITSVertexer;

class AliITSMeanVertexer : public TObject {

 public:
  // default constructor
    AliITSMeanVertexer(Bool_t mode = kTRUE);   
    virtual ~AliITSMeanVertexer();
 
    Bool_t Init();
 
    void   SetFilterOnContributors(Int_t nc=1){fFilterOnContributors = nc;}
    Bool_t Reconstruct(AliRawReader *rawReader);
    void SetCutOnErrX(Double_t cut=3.){fErrXCut = cut;}
    void SetCutOnR(Double_t cut=0.2){fRCut = cut;}
    void SetCutOnCls(UInt_t vmin=50, UInt_t vmax=4000){fLowSPD0=vmin; fHighSPD0=vmax;}
    void   WriteVertices(const char *filename);

    const TH2F*GetVertexXY() const { return fVertexXY; }
    const TH1F*GetVertexZ()  const { return fVertexZ;  }

    // Methods used for debug purposes
    Int_t GetArrayEntries() const {return fVertArray.GetEntries();}
    const AliESDVertex* GetElementAt(Int_t i) const {return (AliESDVertex*)fVertArray[i];}
    UInt_t GetSPD0cls(Int_t i) const {return fClu0[i];}
    Bool_t IsGoodVertex(Int_t i) const {return (fAccEvents.TestBitNumber(i)); }
 
 private:
    // copy constructor (NO copy allowed)
    AliITSMeanVertexer(const AliITSMeanVertexer& vtxr);
    // assignment operator (NO assignment allowed)
    AliITSMeanVertexer& operator=(const AliITSMeanVertexer& /* vtxr */);

    Bool_t Filter(AliESDVertex *vert,UInt_t mult);
    void   AddToMean(AliESDVertex *vert);
    Bool_t ComputeMean(Bool_t killOutliers);
    void Reset();
    void ResetArray(){fAccEvents.ResetAllBits(kTRUE); fVertArray.Clear();
      fIndex=0; for(Int_t i=0;i<fgkMaxNumOfEvents;i++)fClu0[i]=0; }

    static const Int_t fgkMaxNumOfEvents;   // max. number of events 
    AliITSDetTypeRec *fDetTypeRec; //! ITS reco class
    TH2F    *fVertexXY;        //! histogram with transverse vertex distribution (vertex diamond)
    TH1F    *fVertexZ;         //! histogram with longitudinal vertex distribution

    Double_t fWeighPosSum[3];    //! weighted average position sum (transient)
    Double_t fWeighSigSum[3];    //! weighted average position sum (transient)
    Double_t fAverPosSum[3];     //! average position sum (transient)
    Double_t fAverPosSqSum[3][3];//! average square position sum for covariance (transient)
    Double_t fWeighPos[3];       //! weighted average position
    Double_t fWeighSig[3];       //! errors on weighted average position
    Double_t fAverPos[3];        //! average position
    Double_t fAverPosSq[3][3];   //! average square position for covariance
    Int_t fNoEventsContr;        //! number of events used for mean vertex
    Double_t fTotContributors;      //! Integrated number of contributors
    Double_t fAverContributors;  //! Average number of contributors
    Int_t fFilterOnContributors; //! Numb. of contrib must be > fFilter...
    Bool_t fMode;                //! kTRUE for Vertexer3D; 
                                 //! kFALSE for VertexerTapan
    AliITSVertexer* fVertexer;   //! vertexer

    TBits fAccEvents;           //! bit string: 1 for good events 0 for bad ones
    TClonesArray fVertArray;     //! array of vertices to be averaged
    UInt_t *fClu0;              //! number of clusters on SPD inner layer
    Int_t fIndex;               //! current index on the arrays
    Double_t fErrXCut;          //! cut on error on X (error*1000<fErrXCut)
    Double_t fRCut;             //| cut on distance from first estimate (mm)
    UInt_t fLowSPD0;            //! low SPD0 cls value to accept event
    UInt_t fHighSPD0;           //! high SPD0 cls value to accept event
    TH1F *fMultH;               //! debug hist: mult. on SPD0 before Filter
    TH1F *fErrXH;               //! debug hist: error on X before Filter
    TH1F *fMultHa;              //! debug hist: mult. on SPD0 after Filter
    TH1F *fErrXHa;              //! debug hist: error on X after Filter
    TH1F *fDistH;               //! debug hist: distance from peak

    ClassDef(AliITSMeanVertexer,0)
};

#endif
