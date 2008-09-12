#ifndef ALIITSVERTEXER3DTAPAN_H
#define ALIITSVERTEXER3DTAPAN_H

/*  See cxx source for full Copyright notice */


//-------------------------------------------------------------------------
//                          Class AliITSVertexer3DTapan
//   This is a class for the 3d vertex finding
//
//    Origin: Tapan Nayak, VECC-CERN, Tapan.Nayak@cern.ch
//-------------------------------------------------------------------------



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                        AliITSVertexer3DTapan                         //
//                                                                      //
//           Implementation of a 3D vertex finder based on              //
//           SPD clusters.                                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TArrayD.h>
#include "AliITSVertexer.h"

class TTree;
class AliESDVertex;

class AliITSVertexer3DTapan : public AliITSVertexer {
public:
  AliITSVertexer3DTapan(Int_t n=10000):
    fX1(n),fY1(n),fZ1(n), fPhi1(n), ficlu1(0),
    fX2(n),fY2(n),fZ2(n), fPhi2(n), ficlu2(0) {;}
  virtual ~AliITSVertexer3DTapan(){}
  virtual AliESDVertex *FindVertexForCurrentEvent(TTree *cTree);
  virtual void PrintStatus() const {}

protected:
  void LoadClusters(TTree *cTree);
  void CalculatePhi(Float_t fx, Float_t fy, Float_t & phi);
  void CalculateVertex3d1(Double_t pos[3], Float_t cuts[3], Int_t &ncontr);
  void CalculateVertex3d2(Double_t pos[3], Float_t cuts[3], Int_t &ncontr, Double_t sigpos[3]);
  void CalculateLine(Double_t P1[4], Double_t P2[4], Double_t P3[4], Double_t P4[4], Double_t Pa[3], Double_t Pb[3]) const;

private:
  AliITSVertexer3DTapan(AliITSVertexer3DTapan &);
  AliITSVertexer3DTapan& operator=(const AliITSVertexer3DTapan &);

  TArrayD fX1;     // X position of cluster on layer 1 of ITS
  TArrayD fY1;     // Y position of cluster on layer 1 of ITS
  TArrayD fZ1;     // Z position of cluster on layer 1 of ITS
  TArrayD fPhi1;   // Phi position of cluster on layer 1 of ITS
  Int_t   ficlu1;   // Number of clusters on layer 1 of ITS
   
  TArrayD fX2;      // X position of cluster on layer 2 of ITS
  TArrayD fY2;      // Y position of cluster on layer 2 of ITS
  TArrayD fZ2;      // Z position of cluster on layer 2 of ITS
  TArrayD fPhi2;    // Phi position of cluster on layer 2 of ITS
  Int_t   ficlu2;    // Number of clusters on layer 2 of ITS
   
  ClassDef(AliITSVertexer3DTapan,3);
};

#endif
