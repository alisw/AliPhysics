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

#include <TF1.h>

#include <AliAnalysisUtils.h>
#include <AliAnalysisCuts.h>

class AliTriggerAnalysis;
class AliESDVertex;
class AliAODVertex;


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
  void SetCentralityRange(Double_t min, Double_t max) { fCentMin=min; fCentMax=max; }
  void SetCutOnV0MultipicityNTrks(TF1* parMean, TF1* parSigma, Double_t cutSigma=3.) { fparMean=parMean; fparSigma=parSigma; fcutSigma=cutSigma; }
  void SetCutOnNVtxContributorsGloablTPC(TF1* parMin, TF1* parMax) { fparMinVtxContributors=parMin; fparMaxVtxContributors=parMax; }
  void SetRequire2013vertexandevent(Bool_t req13 = kTRUE) {fRequire13sel = req13; }
  //
  //Analysis cuts interface
  //
  virtual Bool_t IsSelected(TObject* event);
  Bool_t IsSelectedESD(TObject* event);
  Bool_t IsSelectedAOD(TObject* event);
  virtual Bool_t IsSelected(TList*   /* list */ ) {return kFALSE;}

  void Print(const Option_t* option = "") const;

private:
  static const char* fgkVtxNames[AliDielectronEventCuts::kVtxTracksOrSPD+1];  //vertex names
  Double_t fVtxZmin;                // minimum z vertex position
  Double_t fVtxZmax;                // maximum z vertex position
  Bool_t   fRequireVtx;             // require a vertex
  Int_t    fMinVtxContributors;     // min number of vertex contributors
  Bool_t   fMultITSTPC;             // if to cut on the ITS TPC multiplicity correlation (Pb-Pb)
  Double_t fCentMin;                // minimum multiplity percentile
  Double_t fCentMax;                // maximum multiplity percentile
  EVtxType fVtxType;                // vertex type
  Bool_t fRequire13sel;             //bit to select event and vertex selection proposed for 2013 in 
                                    //https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PAVertexSelectionStudies
  AliAnalysisUtils fUtils;          //data member to use utility class for event and vertex selection in 2013
  

  UChar_t fRequireV0and;             // use V0and triggered events only

  AliTriggerAnalysis *fTriggerAnalysis; //! trigger analysis class
  const AliESDVertex *fkVertex;         //! current vertex
  const AliAODVertex *fkVertexAOD;      //! current vertex AOD

  TF1*     fparMean;                // parametrization of the mean values
  TF1*     fparSigma;               // parametrization of the sigmas
  Double_t fcutSigma;               // number of absolut sigmas inclusion
  TF1*     fparMinVtxContributors;  // parametrization of #vtx contributors global vs TPC (lower limit) 
  TF1*     fparMaxVtxContributors;  // parametrization of #vtx contributors global vs TPC (upper limit) 
  AliDielectronEventCuts(const AliDielectronEventCuts &c);
  AliDielectronEventCuts &operator=(const AliDielectronEventCuts &c);

  
  ClassDef(AliDielectronEventCuts,2)         // Dielectron EventCuts
};


#endif
