//
// Class AliRsnHistoDef
//
// Definition for a histogram type.
// Since one could do an analysis which is not an invariant mass
// the histogram definition should be more flexible, and it is stored
// separately in a new class.
// This class considers the possibility of a 1D or 2D histograms
// with its related binning, and can create a new histo from his definitions
//

#ifndef AliRsnHistoDef_H
#define AliRsnHistoDef_H

class TH1;
class TH1D;
class TH2D;

class AliRsnHistoDef : public TObject
{
public:

    AliRsnHistoDef();
    AliRsnHistoDef(Int_t n, Double_t min, Double_t max);
    AliRsnHistoDef(Int_t nX, Double_t minX, Double_t maxX, Int_t nY, Double_t minY, Double_t maxY);
    virtual ~AliRsnHistoDef() { }

    Int_t    GetNDimensions() const {return fNDim;}
    Int_t    GetNBins(Int_t axis=0) const {return (axis==0)?fNBins[0]:fNBins[1];}
    Double_t GetMin(Int_t axis=0) const {return (axis==0)?fMin[0]:fMin[1];}
    Double_t GetMax(Int_t axis=0) const {return (axis==0)?fMax[0]:fMax[1];}

    void     SetBins(Int_t n, Double_t min, Double_t max);
    void     SetBins(Int_t nX, Double_t minX, Double_t maxX, Int_t nY, Double_t minY, Double_t maxY);

    TH1D*    Create1DHistogram(const char *name, const char *title);
    TH2D*    Create2DHistogram(const char *name, const char *title);

private:

    void  CheckEdges();

    Int_t        fNDim;       // number of dimensions
    Int_t        fNBins[2];   // number of bins
    Double_t     fMin[2];     // lower edge
    Double_t     fMax[2];     // upper edge
    
    // ROOT dictionary
    ClassDef(AliRsnHistoDef, 1)
};

#endif
