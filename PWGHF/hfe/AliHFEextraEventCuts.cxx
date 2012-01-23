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
// the requirements on the number of charged tracks and on 
// the vertex position and resolution are implemented
// The argument of IsSelected member function (passed object) is cast into 
// an AliESDEvent. In the future may be modified to use AliVEvent interface
// and include more cut variables.
// The class derives from AliCFCutBase
// Author:S.Arcelli Silvia.Arcelli@cern.ch
//
//
#include "TH1F.h"
#include "TList.h"
#include "TBits.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliHFEextraEventCuts.h"
ClassImp(AliHFEextraEventCuts) 
//____________________________________________________________________
AliHFEextraEventCuts::AliHFEextraEventCuts() : 
  AliCFCutBase(),
  fRequireVtxCuts(kFALSE),
  fVtxZMax(1.e99),
  fVtxZMin(-1.e99),
  fVtxNCtrbMin(0),
  fVtxMixed(0),
  fBitMap(0x0)
{
  //
  //ctor
  //
}

//____________________________________________________________________
AliHFEextraEventCuts::AliHFEextraEventCuts(Char_t* name, Char_t* title) : 
  AliCFCutBase(name,title),
  fRequireVtxCuts(kFALSE),
  fVtxZMax(1.e99),
  fVtxZMin(-1.e99),
  fVtxNCtrbMin(0),
  fVtxMixed(0),
  fBitMap(0x0)
 {
  //
  //ctor
  //
  fBitMap=new TBits(0);
  Initialise();
 }

//____________________________________________________________________
AliHFEextraEventCuts::AliHFEextraEventCuts(const AliHFEextraEventCuts& c) : 
  AliCFCutBase(c),
  fRequireVtxCuts(c.fRequireVtxCuts),
  fVtxZMax(c.fVtxZMax),
  fVtxZMin(c.fVtxZMin),
  fVtxNCtrbMin(c.fVtxNCtrbMin),
  fVtxMixed(c.fVtxMixed),
  fBitMap(c.fBitMap)
{
  //
  //copy constructor
  //
  for (Int_t i=0; i<c.kNCuts; i++){
    for (Int_t j=0; j<c.kNStepQA; j++){
      if(c.fhQA[i][j]) fhQA[i][j] = (TH1F*)c.fhQA[i][j]->Clone();
    }
  }

}

//____________________________________________________________________
AliHFEextraEventCuts::~AliHFEextraEventCuts() {
  //
  //dtor
  //

  for (Int_t i=0; i<kNCuts; i++){
    for (Int_t j=0; j<kNStepQA; j++){
      if(fhQA[i][j]) delete fhQA[i][j];
    }
  }

  if(fBitMap)delete fBitMap;

}
//_____________________________________________________________________________
void AliHFEextraEventCuts::Initialise()
{

  //
  //initialization
  //

  //
  // sets pointers to histos to zero
  //

  for(Int_t i=0; i<kNCuts; i++){
    for(Int_t j =0; j<kNStepQA; j++){
      fhQA[i][j]=0x0;
    }
  }
}

//____________________________________________________________________
AliHFEextraEventCuts& AliHFEextraEventCuts::operator=(const AliHFEextraEventCuts& c)
{
  //
  // Assignment operator
  //
  if (this != &c) {
    AliCFCutBase::operator=(c) ;
    fRequireVtxCuts=c.fRequireVtxCuts;
    fVtxZMax=c.fVtxZMax;
    fVtxZMin=c.fVtxZMin;
    fVtxNCtrbMin=c.fVtxNCtrbMin;
    fVtxMixed=c.fVtxMixed;
    fBitMap=c.fBitMap;
  }

  for (Int_t i=0; i<c.kNCuts; i++){
    for (Int_t j=0; j<c.kNStepQA; j++){
      if(c.fhQA[i][j]) fhQA[i][j] = (TH1F*)c.fhQA[i][j]->Clone();
    }
  }


  return *this ;
}

//____________________________________________________________________
Bool_t AliHFEextraEventCuts::IsSelected(TObject* obj) {
  //
  //Check if the requested cuts are passed
  //


  SelectionBitMap(obj);

  if (fIsQAOn) FillHistograms(obj,0);
  Bool_t isSelected = kTRUE;

  for (UInt_t icut=0; icut<fBitMap->GetNbits();icut++)
	if(!fBitMap->TestBitNumber(icut)) isSelected = kFALSE;

  if (!isSelected) return kFALSE ;
  if (fIsQAOn) FillHistograms(obj,1);
  return kTRUE;

}
//____________________________________________________________________
void AliHFEextraEventCuts::SelectionBitMap(TObject* obj) {
  //
  //cut on the number of charged tracks and on the event vertex.
  //so far specific to AliESDEvents
  //

  //Check if the requested cuts are passed and return a bitmap
  for(Int_t j=0;j<kNCuts;j++)fBitMap->SetBitNumber(j,kFALSE);
  AliESDEvent* esd = dynamic_cast<AliESDEvent *>(obj);
  if ( !esd ) return;

  //now start checking the cuts,
  //first assume the event will be accepted: 
  for(Int_t j=0;j<kNCuts;j++)fBitMap->SetBitNumber(j,kTRUE);


  
  if(fRequireVtxCuts){
    const AliESDVertex* vtxESD = 0x0;
    if      (fVtxMixed) {
      vtxESD = esd->GetPrimaryVertexTracks() ;
      if((!vtxESD) || (vtxESD->GetNContributors() <= 0)) {
	vtxESD = esd->GetPrimaryVertexSPD() ;
	if((!vtxESD) || (vtxESD->GetNContributors() <= 0)) {
	  for(Int_t j=1;j<kNCuts;j++)fBitMap->SetBitNumber(j,kFALSE); 
	  AliWarning("Cannot get vertex, skipping event");
	  return;
	}
      }
    }
    else {   
      vtxESD = esd->GetPrimaryVertexTracks() ;
    }

    if(!vtxESD){
      for(Int_t j=1;j<kNCuts;j++)fBitMap->SetBitNumber(j,kFALSE); 
      AliWarning("Cannot get vertex, skipping event");
      return;
    }
    
    // Pick up the position and uncertainties
    vtxESD = esd->GetPrimaryVertex();
    if(!vtxESD){
      for(Int_t j=1;j<kNCuts;j++)fBitMap->SetBitNumber(j,kFALSE); 
      AliWarning("Cannot get vertex, skipping event");
      return;
    }
    Double_t vtxPos[3];
    vtxPos[0] = vtxESD->GetXv();
    vtxPos[1] = vtxESD->GetYv();
    vtxPos[2] = vtxESD->GetZv();
    
    Int_t nCtrb = vtxESD->GetNContributors();

    // Apply the cut
    
    if (vtxPos[2]>fVtxZMax || vtxPos[2]<fVtxZMin)
      fBitMap->SetBitNumber(0,kFALSE); 
    if (nCtrb<fVtxNCtrbMin)
      fBitMap->SetBitNumber(1,kFALSE);

  }  
  return;
}

//_____________________________________________________________________________
void AliHFEextraEventCuts::FillHistograms(TObject* obj, Bool_t b)
{
  //
  // fill the QA histograms
  //

  if(!fIsQAOn) return;
  // cast TObject into VParticle
  AliESDEvent* esd = dynamic_cast<AliESDEvent *>(obj);
  if (!esd ) return  ;

  // index = 0: fill histograms before cuts
  // index = 1: fill histograms after cuts
  Int_t index = -1;
  index = ((b) ? 1 : 0);


  //look at vertex parameters:
  const AliESDVertex* vtxESD = 0x0;
  if      (fVtxMixed) {
    vtxESD = esd->GetPrimaryVertexTracks() ;
    if((!vtxESD) || (vtxESD->GetNContributors() <= 0)) {
      vtxESD = esd->GetPrimaryVertexSPD() ;
      if((!vtxESD) || (vtxESD->GetNContributors() <= 0)) {
	return;
      }
    }
  }
  else {   
    vtxESD = esd->GetPrimaryVertexTracks() ;
  }
  if(!vtxESD)return;
  // vertex position and uncertainties
  fhQA[kVtxPosZ] [index]->Fill(vtxESD->GetZv());
  fhQA[kVtxNCtrb][index]->Fill(vtxESD->GetNContributors());
  
}

//____________________________________________________________________
void AliHFEextraEventCuts::SetHistogramBins(Int_t index, Int_t nbins, Double_t *bins)
{
  //
  //setting x-axis bin limits of QA histogram fhQA[index] 
  // 

  for(Int_t i=0;i<kNStepQA;i++){
    if(!fhQA[index][i]){AliWarning("non-existing histogram!");
    return;
    }
    fhQA[index][i]->GetXaxis()->Set(nbins,bins);
  }
}
//____________________________________________________________________
void AliHFEextraEventCuts::SetHistogramBins(Int_t index, Int_t nbins, Double_t xmin, Double_t xmax)
{
  //
  //setting x-axis bins and range of QA histogram fhQA[index] 
  // 

  for(Int_t i=0;i<kNStepQA;i++){
    if(!fhQA[index][i]){AliWarning("non-existing histogram!");
    return;
    }
    fhQA[index][i]->GetXaxis()->Set(nbins,xmin,xmax);
  }
}

//_____________________________________________________________________________
 void AliHFEextraEventCuts::DefineHistograms() {
  //
  // histograms for cut variables
  //
  Int_t color = 2;

  if(!fIsQAOn) {
    AliInfo(Form("No QA histos requested, Please first set the QA flag on!"));
    return;
  }  
  
  // book QA histograms

  Char_t str[5];
  for (Int_t i=0; i<kNStepQA; i++) {
    if (i==0) snprintf(str,5," ");
    else snprintf(str,5,"_cut");

    fhQA[kVtxPosZ][i]	= new  TH1F(Form("%s_Vtx_Pos_Z%s",GetName(),str),		"",200,-50.,50.);
    fhQA[kVtxNCtrb][i]	= new  TH1F(Form("%s_Vtx_N_Ctrb%s",GetName(),str),		"",1000,0.,1000.);
 
    fhQA[kVtxPosZ][i]	->SetXTitle("Vertex Position Z (cm)");
    fhQA[kVtxNCtrb][i]	->SetXTitle("Number of contributors");
  }

  for(Int_t i=0; i<kNCuts; i++) fhQA[i][1]->SetLineColor(color);

}

//_____________________________________________________________________________
void AliHFEextraEventCuts::AddQAHistograms(TList *qaList) {
  //
  // saves the histograms in a TList
  //

  DefineHistograms();

  for (Int_t j=0; j<kNStepQA; j++) {
    for(Int_t i=0; i<kNCuts; i++)
	qaList->Add(fhQA[i][j]);
  }
}
