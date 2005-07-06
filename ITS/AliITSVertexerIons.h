#ifndef ALIITSVERTEXERIONS_H
#define ALIITSVERTEXERIONS_H

#include <AliITSVertexer.h>

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


class TH1F;

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

  Int_t fNpThreshold;      // minimum number of rec points for vertexing
  Double_t fMaxDeltaPhi;   // Maximum phi difference for rec points correlation
  Double_t fMaxDeltaZ;     // Maximum z difference for rec points correlation
  AliITSVertexerIons(const AliITSVertexerIons &source); // copy constructor (NO copy allowed: the constructor is protected to avoid misuse)   
  AliITSVertexerIons& operator=(const AliITSVertexerIons &source); // assignment operator (NO assignment allowed)

  ClassDef(AliITSVertexerIons,4);
};

#endif
