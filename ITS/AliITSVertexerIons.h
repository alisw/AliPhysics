#ifndef ALIITSVERTEXERIONS_H
#define ALIITSVERTEXERIONS_H

#include <AliITSVertexer.h>

//////////////////////////////////////////////////////////////////////
// AliITSVertexerIons is a class for full 3D primary vertex         //
// finding optimized for Ion-Ion interactions                       //
//                                                                  // 
//                                                                  //
//                                                                  //
//                                                                  //
// Written by Giuseppe Lo Re and Francesco Riggi                    //
// Giuseppe.Lore@ct.infn.it                                         //
//                                                                  //
// Franco.Riggi@ct.infn.it                                          //
//                                                                  //
// Release date: May 2001                                           //
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
  virtual Double_t PhiFunc(Float_t p[]);
  virtual void PrintStatus() const;
  Int_t GetNpThreshold() const {return fNpThreshold;}
  void SetNpThreshold(Int_t t = 4000){fNpThreshold = t;}

 protected:
  AliITS *fITS;            //! pointer to the AliITS object
  Int_t fNpThreshold;      // minimum number of rec points for vertexing

  ClassDef(AliITSVertexerIons,2);
};

#endif
