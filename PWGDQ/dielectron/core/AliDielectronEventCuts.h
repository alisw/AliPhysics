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
#include <TVectorD.h>
#include <TBits.h>

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

  void SetRunRejection(const TVectorD * vec)    { fRun.Use(vec->GetNrows(),vec->GetMatrixArray()); }
  void SetVertexType(EVtxType type)             { fVtxType=type;                }
  void SetVertexZ(Double_t zmin, Double_t zmax) { fVtxZmin=zmin; fVtxZmax=zmax; }
  void SetRequireVertex(Bool_t req=kTRUE)       { fRequireVtx=req;              }
  void SetMaxVertexDifference(Float_t diff = 0.2)    { fVtxDiff = diff;              }
  void SetRequireV0and(UChar_t type=1)          { fRequireV0and=type;           }
  void SetRequireCompleteDAQ(Bool_t req =kTRUE) { f2015IsIncompleteDAQ=req;     }
  void SetMinVtxContributors(Int_t min=1)       { fMinVtxContributors=min;      }
  void SetCutOnMultipicityITSTPC(Bool_t mult=kTRUE) { fMultITSTPC=mult;         }
  void SetCentralityRange(Double_t min, Double_t max, Bool_t userun2=kFALSE) { fCentMin=min; fCentMax=max;fRun2=userun2;}
  void SetCutOnV0MultipicityNTrks(TF1* parMean, TF1* parSigma, Double_t cutSigma=3.) { fparMean=parMean; fparSigma=parSigma; fcutSigma=cutSigma; }
  void SetCutOnNVtxContributorsGloablTPC(TF1* parMin, TF1* parMax) { fparMinVtxContributors=parMin; fparMaxVtxContributors=parMax; }
  void SetRequire2013vertexandevent(Bool_t req13 = kTRUE) {fRequire13sel = req13; }
  void SetMinCorrCutFunction(TF1 *fun, UInt_t varx, UInt_t vary=0);
  void SetMaxCorrCutFunction(TF1 *fun, UInt_t varx, UInt_t vary=0);

  //
  //Analysis cuts interface
  //
  virtual Bool_t IsSelected(TObject* event);
  Bool_t IsSelectedESD(TObject* event);
  Bool_t IsSelectedAOD(TObject* event);
  virtual Bool_t IsSelected(TList*   /* list */ ) {return kFALSE;}


  enum ETRDTriggerClass{ kSE, kQU, kSEorQU, kSEandQU};
  static Bool_t IsTRDTriggerFired( const AliVEvent* event, const ETRDTriggerClass triggerClass, Bool_t &trackMatched, Int_t &bin );

  void Print(const Option_t* option = "") const;

private:
  static const char* fgkVtxNames[AliDielectronEventCuts::kVtxTracksOrSPD+1];  //vertex names

  TBits     *fUsedVars;             // list of used variables
  TVectorD fRun;                    // run rejection vector
  Double_t fVtxZmin;                // minimum z vertex position
  Double_t fVtxZmax;                // maximum z vertex position
  Bool_t   fRequireVtx;             // require a vertex
  Int_t    fMinVtxContributors;     // min number of vertex contributors
  Bool_t   fMultITSTPC;             // if to cut on the ITS TPC multiplicity correlation (Pb-Pb)
  Double_t fCentMin;                // minimum multiplicity percentile
  Double_t fCentMax;                // maximum multiplicity percentile
  Bool_t fRun2;                     //using run2 centrality
  EVtxType fVtxType;                // vertex type
  Bool_t fRequire13sel;             //bit to select event and vertex selection proposed for 2013 in
                                    //https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PAVertexSelectionStudies
  Bool_t f2015IsIncompleteDAQ;	    // Needed for Run2-2015 data analysis, where a problem with incomplete events from the daq-side exists -- kFALSE will reject incomplete events
  Float_t fVtxDiff;                 // cut values for a cut on the difference between the position of the vtx from SPD and tracks. Used in LHC15o std. 0.2 plus additional cut on both vtx resolutions
  AliAnalysisUtils fUtils;          //data member to use utility class for event and vertex selection in 2013


  UChar_t fRequireV0and;             // use V0and triggered events only

  AliTriggerAnalysis *fTriggerAnalysis; //! trigger analysis class
  const AliESDVertex *fkVertex;         //! current vertex
  const AliAODVertex *fkVertexAOD;      //! current vertex AOD

  TH1D* fCorrCutMin[5];       //parametrization of lower limit correlation cut
  TH1D* fCorrCutMax[5];       //parametrization of upper limit correlation cut

  TF1*     fparMean;                // parametrization of the mean values
  TF1*     fparSigma;               // parametrization of the sigmas
  Double_t fcutSigma;               // number of absolut sigmas inclusion
  TF1*     fparMinVtxContributors;  // parametrization of #vtx contributors global vs TPC (lower limit)
  TF1*     fparMaxVtxContributors;  // parametrization of #vtx contributors global vs TPC (upper limit)
  AliDielectronEventCuts(const AliDielectronEventCuts &c);
  AliDielectronEventCuts &operator=(const AliDielectronEventCuts &c);


  ClassDef(AliDielectronEventCuts,4)         // Dielectron EventCuts
};


#endif
