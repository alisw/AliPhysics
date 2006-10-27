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
    virtual ~AliStrLine(); // destructor
    void PrintStatus() const;
    void SetP0(Double_t *point) {for(Int_t i=0;i<3;i++)fP0[i]=point[i];}
    void SetCd(Double_t *cd) {for(Int_t i=0;i<3;i++)fCd[i]=cd[i];}
    void SetDebug(Int_t dbfl = 0){fDebug = dbfl; }  
    void GetP0(Double_t *point) const {for(Int_t i=0;i<3;i++)point[i]=fP0[i];}
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
    Double_t fCd[3];           // direction cosines
    Double_t fTpar;            //! parameter 
    Int_t   fDebug;           //! debug flag - verbose printing if >0
 private:
    void SetPar(Double_t par){fTpar = par;}

  ClassDef(AliStrLine,1);
};

#endif
