#ifndef ALIDIELECTRONEVENTCUTS_H
#define ALIDIELECTRONEVENTCUTS_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#############################################################
//#                                                           # 
//#         Class AliDielectronEventCuts                     #
//#                                                           #
//#  Authors:                                                 #
//#   Anton     Andronic, GSI / A.Andronic@gsi.de             #
//#   Ionut C.  Arsene,   GSI / I.C.Arsene@gsi.de             #
//#   Julian    Book,     Uni Ffm / Julian.Book@cern.ch       #
//#   Frederick Kramer,   Uni Ffm, / Frederick.Kramer@cern.ch #
//#   Magnus    Mager,    CERN / Magnus.Mager@cern.ch         #
//#   WooJin J. Park,     GSI / W.J.Park@gsi.de               #
//#   Jens      Wiechula, Uni HD / Jens.Wiechula@cern.ch      #
//#                                                           #
//#############################################################

#include <AliAnalysisCuts.h>

class AliTriggerAnalysis;
class AliESDVertex;

class AliDielectronEventCuts : public AliAnalysisCuts {
public:
  enum EVtxType { kVtxTracks=0, kVtxSPD, kVtxTPC, kVtxAny, kVtxTracksOrSPD };
  
  AliDielectronEventCuts();
  AliDielectronEventCuts(const char*name, const char* title);

  virtual ~AliDielectronEventCuts();


  void SetVertexType(EVtxType type)             { fVtxType=type;                }
  void SetVertexZ(Double_t zmin, Double_t zmax) { fVtxZmin=zmin; fVtxZmax=zmax; }
  void SetRequireVertex(Bool_t req=kTRUE)       { fRequireVtx=req;              }
  void SetRequireV0and(UChar_t type=1)          { fRequireV0and=type;           }
  void SetMinVtxContributors(Int_t min=1)       { fMinVtxContributors=min;      }
  void SetCutOnMultipicityITSTPC(Bool_t mult=kTRUE) { fMultITSTPC=mult;         }
  //
  //Analysis cuts interface
  //
  virtual Bool_t IsSelected(TObject* event);
  virtual Bool_t IsSelected(TList*   /* list */ ) {return kFALSE;}
  

private:

  Double_t fVtxZmin;                // minimum z vertex position
  Double_t fVtxZmax;                // maximum z vertex position
  Bool_t   fRequireVtx;             // require a vertex
  Int_t    fMinVtxContributors;     // min number of vertex contributors
  Bool_t   fMultITSTPC;             // if to cut on the ITS TPC multiplicity correlation (Pb-Pb)
  EVtxType fVtxType;                // vertex type

  UChar_t fRequireV0and;             // use V0and triggered events only

  AliTriggerAnalysis *fTriggerAnalysis; //! trigger analysis class
  const AliESDVertex *fkVertex;         //! current vertex

  AliDielectronEventCuts(const AliDielectronEventCuts &c);
  AliDielectronEventCuts &operator=(const AliDielectronEventCuts &c);

  
  ClassDef(AliDielectronEventCuts,1)         // Dielectron EventCuts
};



#endif
