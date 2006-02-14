#ifndef ALIITSSTRLINE_H
#define ALIITSSTRLINE_H
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////
//                                                               //
// A straight line is coded as a point (3 Double_t) and           //
// 3 direction cosines                                           //
//                                                               //
///////////////////////////////////////////////////////////////////


class AliITSStrLine : public TObject {

 public:
    AliITSStrLine();        // default constructor
    AliITSStrLine(Double_t *point, Double_t *cd);  // standard constructor
    virtual ~AliITSStrLine(); // destructor
    void PrintStatus() const;
    void SetP0(Double_t *point) {for(Int_t i=0;i<3;i++)fP0[i]=point[i];}
    void SetCd(Double_t *cd) {for(Int_t i=0;i<3;i++)fCd[i]=cd[i];}
    void SetDebug(Int_t dbfl = 0){fDebug = dbfl; }  
    void GetP0(Double_t *point) const {for(Int_t i=0;i<3;i++)point[i]=fP0[i];}
    void GetCd(Double_t *cd) const {for(Int_t i=0;i<3;i++)cd[i]=fCd[i];}
    void GetCurrentPoint(Double_t *point) const;
    Int_t IsParallelTo(AliITSStrLine *line) const;
    Int_t Crossrphi(AliITSStrLine *line);
    Int_t CrossPoints(AliITSStrLine *line, Double_t *point1, Double_t *point2);
    Int_t Cross(AliITSStrLine *line, Double_t *point);
    Double_t GetDCA(AliITSStrLine *line);
 protected:
    Double_t fP0[3];           // given point
    Double_t fCd[3];           // direction cosines
    Double_t fTpar;            //! parameter 
    Int_t   fDebug;           //! debug flag - verbose printing if >0
 private:
    void SetPar(Double_t par){fTpar = par;}

  ClassDef(AliITSStrLine,1);
};

#endif
