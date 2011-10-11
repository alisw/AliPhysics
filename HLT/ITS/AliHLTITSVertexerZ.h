//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTITSVERTEXERZ_H
#define ALIHLTITSVERTEXERZ_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                          High Level Trigger ITS vertexer
//       This class is a fast version of the off-line AliITSVertexerZ.
//       The two main differences with respect to the off-line vertexer
//       are the splitting of the clusters in phi bins and the filling
//       of local arrays instead of root histograms.
//      
//           Origin: Cvetan Cheshkov, CERN, Cvetan.Cheshkov@cern.ch 
//-------------------------------------------------------------------------

#include "AliITSVertexerZ.h"

class TString;
class TTree;
class AliITSgeom;

//-------------------------------------------------------------------------
class AliHLTITSVertexerZ : public AliITSVertexerZ {
public:
  AliHLTITSVertexerZ();
  AliHLTITSVertexerZ(Float_t x0, Float_t y0);
  virtual ~AliHLTITSVertexerZ();

  AliESDVertex* FindVertexForCurrentEvent(AliITSgeom* /* geom */,TTree *tR);

  void SetBinWidthFine(Float_t bw=0.0005){fStepFine = bw;}

 private:
  AliHLTITSVertexerZ(const AliHLTITSVertexerZ &vtxr);
  AliHLTITSVertexerZ& operator=(const AliHLTITSVertexerZ&  vtxr );

  TH1F *fZCombf;           //! histogram with fine z distribution
  Float_t fStepFine;       // bin width for fZCombf histogram

  ClassDef(AliHLTITSVertexerZ,3)   //HLT ITS vertexer
};

typedef AliHLTITSVertexerZ AliL3ITSVertexerZ; // for backward compatibility

#endif
