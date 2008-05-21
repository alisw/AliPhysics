#ifndef ALISTRLINE_H
#define ALISTRLINE_H
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////
//                                                               //
// A straight line is coded as a point (3 Double_t) and           //
// 3 direction cosines                                           //
//                                                               //
///////////////////////////////////////////////////////////////////


class AliStrLine : public TObject {

 public:
    AliStrLine();        // default constructor
    AliStrLine(Double_t *point, Double_t *cd, Bool_t twopoints=kFALSE);  // standard constructor
    AliStrLine(Float_t *pointf, Float_t *cdf, Bool_t twopoints=kFALSE); 
    AliStrLine(Double_t *point, Double_t *sig2point, Double_t *cd, Bool_t twopoints=kFALSE);
    AliStrLine(Float_t *pointf, Float_t *sig2point, Float_t *cdf, Bool_t twopoints=kFALSE); 
    AliStrLine(Double_t *point, Double_t *sig2point, Double_t *wmat, Double_t *cd, Bool_t twopoints=kFALSE);
    AliStrLine(Float_t *pointf, Float_t *sig2point, Float_t *wmat, Float_t *cdf, Bool_t twopoints=kFALSE); 
    AliStrLine(const AliStrLine& source);
    AliStrLine& operator=(const AliStrLine& source);
    virtual ~AliStrLine(); // destructor
    virtual void Clear(Option_t*){if(fWMatrix)delete[] fWMatrix; fWMatrix = 0;}
    void PrintStatus() const;
    void SetP0(const Double_t *point) {for(Int_t i=0;i<3;i++)fP0[i]=point[i];}
    void SetSigma2P0(const Double_t *sigsq) {for(Int_t i=0;i<3;i++)fSigma2P0[i]=sigsq[i];}
    void SetWMatrix(const Double_t *wmat);
    void SetCd(const Double_t *cd) {for(Int_t i=0;i<3;i++)fCd[i]=cd[i];}
    void GetP0(Double_t *point) const {for(Int_t i=0;i<3;i++)point[i]=fP0[i];}
    void GetSigma2P0(Double_t *sigsq) const {for(Int_t i=0;i<3;i++)sigsq[i]=fSigma2P0[i];}
    void GetWMatrix(Double_t *wmat) const;
    void GetCd(Double_t *cd) const {for(Int_t i=0;i<3;i++)cd[i]=fCd[i];}
    void GetCurrentPoint(Double_t *point) const;
    Int_t IsParallelTo(AliStrLine *line) const;
    Int_t Crossrphi(AliStrLine *line);
    Int_t CrossPoints(AliStrLine *line, Double_t *point1, Double_t *point2);
    Int_t Cross(AliStrLine *line, Double_t *point);
    Double_t GetDCA(AliStrLine *line) const;
    Double_t GetDistFromPoint(Double_t *point) const;
 protected:
    void InitDirection(Double_t *point, Double_t *cd);
    void InitTwoPoints(Double_t *pA, Double_t *pB);
    Double_t fP0[3];           // given point
    Double_t fSigma2P0[3];           // errors on coordinates of given point
    Double_t *fWMatrix;           //[6] weighting matrix
    /* fWMatrix is a symmetric matrix internally stored as
       0 --> row = 0, col = 0
       1 --> 0,1
       2 --> 0,2
       3 --> 1,1
       4 --> 1,2
       5 --> 2,2
       The external interface (constructor, getter and setter) is:
       0 --> row = 0, col = 0
       1 --> 0,1
       2 --> 0,2
       3 --> 1,0
       4 --> 1,1
       5 --> 1,2
       6 --> 2,0
       7 --> 2,1
       8 --> 2,2                                                 */
    Double_t fCd[3];           // direction cosines
    Double_t fTpar;            //! parameter 
 private:
    void SetPar(const Double_t par){fTpar = par;}

  ClassDef(AliStrLine,4);
};

#endif
