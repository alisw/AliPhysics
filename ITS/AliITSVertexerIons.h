#ifndef ALIITSVERTEXERIONS_H
#define ALIITSVERTEXERIONS_H

#include <AliITSVertexer.h>
#include <TH1F.h>

//////////////////////////////////////////////////////////////////////
// AliITSVertexerIons  is a class for full 3D primary vertex        //
// finding optimized for Ion-Ion interactions                       //
//                                                                  // 
//                                                                  //
//                                                                  //
//                                                                  //
// Written by Giuseppe Lo Re and Francesco Riggi                    //
// Giuseppe.Lore@ct.infn.it                                         //
// Franco.Riggi@ct.infn.it                                          //
//                                                                  //
// Release date: Mar 2004                                           //
//                                                                  //
//                                                                  //       
//////////////////////////////////////////////////////////////////////

class AliITS;

class AliITSVertexerIons : public AliITSVertexer {

 public:
  AliITSVertexerIons();
  AliITSVertexerIons(TString fn); 
  virtual ~AliITSVertexerIons(); // destructor
  virtual AliESDVertex* FindVertexForCurrentEvent(Int_t event);
  virtual void FindVertices();
  virtual void PhiFunc(Double_t &x,Double_t &y,Double_t &phi);
  virtual void PrintStatus() const;
  Int_t GetNpThreshold() const {return fNpThreshold;}
  void SetNpThreshold(Int_t t = 500){fNpThreshold = t;}
  Double_t GetMaxDeltaPhi() const {return fMaxDeltaPhi;}
  void SetMaxDeltaPhi(Double_t dphi=0.45) {fMaxDeltaPhi=dphi;}
  Double_t GetMaxDeltaZ() const {return fMaxDeltaPhi;}
  void SetMaxDeltaZ(Double_t dz=0.15) {fMaxDeltaZ=dz;}
  Double_t FindMaxAround(Double_t point, TH1F *h, Double_t distance);

 protected:
  AliITS *fITS;            //! pointer to the AliITS object
  Int_t fNpThreshold;      // minimum number of rec points for vertexing
  Double_t fMaxDeltaPhi;
  Double_t fMaxDeltaZ;

  ClassDef(AliITSVertexerIons,3);
};

#endif
