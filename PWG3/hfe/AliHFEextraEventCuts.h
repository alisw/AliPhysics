#ifndef ALIHFEEXTRAEVENTCUTS_H
#define ALIHFEEXTRAEVENTCUTS_H
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
// Cut on the Event at reconstructed level: for the moment 
// just the requirements on the number of charged tracks and on 
// the vertex position and resolution are implemented
// The argument of IsSelected member function (passed object) is cast into 
// an AliESDEvent. In the future may be modified to use AliVEvent interface
// and include more cut variables.
// The class derives from AliCFCutBase
// Author:R.Bailhache


#include "AliCFCutBase.h"
class TH1F;
class TBits;
//_____________________________________________________________________________
class AliHFEextraEventCuts: public AliCFCutBase 
{
 public :
  AliHFEextraEventCuts() ;
  AliHFEextraEventCuts(Char_t* name, Char_t* title) ;
  AliHFEextraEventCuts(const AliHFEextraEventCuts& c) ;
  AliHFEextraEventCuts& operator=(const AliHFEextraEventCuts& c) ;
  ~AliHFEextraEventCuts();
  Bool_t IsSelected(TObject* obj);
  Bool_t IsSelected(TList* /*list*/) {return kTRUE;}

  void SetRequireVtxCuts(Bool_t vtx=kFALSE) {fRequireVtxCuts=vtx;} // cut values setter
  void SetVertexZCut(Double_t zMin=-1.e99, Double_t zMax=1.e99) { fVtxZMin=zMin; fVtxZMax=zMax;} // cut values setter

  void SetVertexNContributors(Int_t min) {fVtxNCtrbMin=min;}
  void SetUseMixedVertex() {fVtxMixed=kTRUE ;} //default is vertex from tracks
 
  Bool_t   GetRequireVtxCuts() const {return fRequireVtxCuts;} // cut value getter
  Double_t GetVertexZMax() const {return fVtxZMax;} // cut values getter
  Double_t GetVertexZMin() const {return fVtxZMin;} // cut values getter
  
  // QA histogram setter
  // please use indices from the enumeration below
  void SetHistogramBins(Int_t index, Int_t nbins, Double_t *bins);
  void SetHistogramBins(Int_t index, Int_t nbins, Double_t xmin, Double_t xmax);
  enum{kVtxPosZ,
       kVtxNCtrb,
       kNCuts,
       kNStepQA=2
  };
  
 protected:
  void SelectionBitMap(TObject* obj);
  void AddQAHistograms(TList *qaList) ;
  void DefineHistograms(); 		// books histograms 
  void Initialise();			// sets everything to 0
  void FillHistograms(TObject* obj, Bool_t b);

  Bool_t fRequireVtxCuts ; //The type of trigger to be checked
  Double_t fVtxZMax ; //Z vertex position, maximum value
  Double_t fVtxZMin ; //Z vertex position, minimum value
  Int_t    fVtxNCtrbMin; //Min number of contributors to vertex
  Bool_t   fVtxMixed;  //Flag for use of mixed vertex (primary vertex with track, if not SPD vertex)
 
  TBits *fBitMap ; //cut mask

  TH1F* fhQA[kNCuts][kNStepQA];		// QA Histograms

  ClassDef(AliHFEextraEventCuts,3);
};

#endif
