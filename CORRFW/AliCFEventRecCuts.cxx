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
#include "AliCFEventRecCuts.h"
ClassImp(AliCFEventRecCuts) 
//____________________________________________________________________
AliCFEventRecCuts::AliCFEventRecCuts() : 
  AliCFCutBase(),
  fNTracksMin(-1),
  fNTracksMax(1000000),
  fRequireVtxCuts(kFALSE),
  fVtxXMax(1.e99),
  fVtxYMax(1.e99),
  fVtxZMax(1.e99),
  fVtxXMin(-1.e99),
  fVtxYMin(-1.e99),
  fVtxZMin(-1.e99),
  fVtxXResMax(1.e99),
  fVtxYResMax(1.e99),
  fVtxZResMax(1.e99),
  fVtxNCtrbMin(0),
  fVtxNCtrbMax((Int_t)1.e9),
  fVtxTPC(0),
  fBitMap(0x0)
{
  //
  //ctor
  //
  fBitMap=new TBits(0);
  Initialise();
}

//____________________________________________________________________
AliCFEventRecCuts::AliCFEventRecCuts(Char_t* name, Char_t* title) : 
  AliCFCutBase(name,title),
  fNTracksMin(-1),
  fNTracksMax(1000000),
  fRequireVtxCuts(kFALSE),
  fVtxXMax(1.e99),
  fVtxYMax(1.e99),
  fVtxZMax(1.e99),
  fVtxXMin(-1.e99),
  fVtxYMin(-1.e99),
  fVtxZMin(-1.e99),
  fVtxXResMax(1.e99),
  fVtxYResMax(1.e99),
  fVtxZResMax(1.e99),
  fVtxNCtrbMin(0),
  fVtxNCtrbMax((Int_t)1.e9),
  fVtxTPC(0),
  fBitMap(0x0)
 {
  //
  //ctor
  //
  fBitMap=new TBits(0);
  Initialise();
 }

//____________________________________________________________________
AliCFEventRecCuts::AliCFEventRecCuts(const AliCFEventRecCuts& c) : 
  AliCFCutBase(c),
  fNTracksMin(c.fNTracksMin),
  fNTracksMax(c.fNTracksMax),
  fRequireVtxCuts(c.fRequireVtxCuts),
  fVtxXMax(c.fVtxXMax),
  fVtxYMax(c.fVtxYMax),
  fVtxZMax(c.fVtxZMax),
  fVtxXMin(c.fVtxXMin),
  fVtxYMin(c.fVtxYMin),
  fVtxZMin(c.fVtxZMin),
  fVtxXResMax(c.fVtxXResMax),
  fVtxYResMax(c.fVtxYResMax),
  fVtxZResMax(c.fVtxZResMax),
  fVtxNCtrbMin(c.fVtxNCtrbMin),
  fVtxNCtrbMax(c.fVtxNCtrbMax),
  fVtxTPC(c.fVtxTPC),
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
AliCFEventRecCuts::~AliCFEventRecCuts() {
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
void AliCFEventRecCuts::Initialise()
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
AliCFEventRecCuts& AliCFEventRecCuts::operator=(const AliCFEventRecCuts& c)
{
  //
  // Assignment operator
  //
  if (this != &c) {
    AliCFCutBase::operator=(c) ;
    fNTracksMin=c.fNTracksMin;
    fNTracksMax=c.fNTracksMax;
    fRequireVtxCuts=c.fRequireVtxCuts;
    fVtxXMax=c.fVtxXMax;
    fVtxYMax=c.fVtxYMax;
    fVtxZMax=c.fVtxZMax;
    fVtxXMin=c.fVtxXMin;
    fVtxYMin=c.fVtxYMin;
    fVtxZMin=c.fVtxZMin;
    fVtxXResMax=c.fVtxXResMax;
    fVtxYResMax=c.fVtxYResMax;
    fVtxZResMax=c.fVtxZResMax;
    fVtxNCtrbMin=c.fVtxNCtrbMin;
    fVtxNCtrbMax=c.fVtxNCtrbMax;
    fVtxTPC=c.fVtxTPC;
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
Bool_t AliCFEventRecCuts::IsSelected(TObject* obj) {
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
void AliCFEventRecCuts::SelectionBitMap(TObject* obj) {
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

  //Number of charged tracks:
  Int_t nTracks = esd->GetNumberOfTracks();
  if(nTracks<fNTracksMin || nTracks>fNTracksMax)
    fBitMap->SetBitNumber(0,kFALSE); 
  
  if(fRequireVtxCuts){
    const AliESDVertex* vtxESD = fVtxTPC ? esd->GetPrimaryVertexTPC() : esd->GetPrimaryVertexSPD() ;
    if(!vtxESD){
      for(Int_t j=1;j<kNCuts;j++)fBitMap->SetBitNumber(j,kFALSE); 
      AliWarning("Cannot get vertex, skipping event");
      return;
    }
    // Require the vertex to have been reconstructed successfully
    if (strcmp(vtxESD->GetName(), "default")==0){
      AliWarning(Form(" No reconstructed vertex found, skipping event"));    
      for(Int_t j=1;j<kNCuts;j++)fBitMap->SetBitNumber(j,kFALSE); 
      return;
    }    
    // Pick up the position and uncertainties
    
    Double_t vtxPos[3];
    vtxPos[0] = vtxESD->GetXv();
    vtxPos[1] = vtxESD->GetYv();
    vtxPos[2] = vtxESD->GetZv();
    
    Double_t vtxRes[3];
    vtxRes[0] = vtxESD->GetXRes();
    vtxRes[1] = vtxESD->GetYRes();
    vtxRes[2] = vtxESD->GetZRes();
 
    Int_t nCtrb = vtxESD->GetNContributors();

    // Apply the cut
    
    if (vtxPos[0]>fVtxXMax || vtxPos[0]<fVtxXMin)
      fBitMap->SetBitNumber(1,kFALSE); 
    if (vtxPos[1]>fVtxYMax || vtxPos[1]<fVtxYMin)
      fBitMap->SetBitNumber(2,kFALSE); 
    if (vtxPos[2]>fVtxZMax || vtxPos[2]<fVtxZMin)
      fBitMap->SetBitNumber(3,kFALSE); 
    if (vtxRes[0]==0 || vtxRes[0]>fVtxXResMax)
      fBitMap->SetBitNumber(4,kFALSE); 
    if (vtxRes[1]==0 || vtxRes[1]>fVtxYResMax)
      fBitMap->SetBitNumber(5,kFALSE); 
    if (vtxRes[2]==0 || vtxRes[2]>fVtxZResMax)
      fBitMap->SetBitNumber(6,kFALSE); 
    if (nCtrb<fVtxNCtrbMin || nCtrb>fVtxNCtrbMax)
      fBitMap->SetBitNumber(7,kFALSE);

  }  
  return;
}

//_____________________________________________________________________________
void AliCFEventRecCuts::FillHistograms(TObject* obj, Bool_t b)
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


  //number of charged tracks:
  Int_t nTracks = esd->GetNumberOfTracks();
  fhQA[kNTracks][index]->Fill(nTracks);

  //look at vertex parameters:
  const AliESDVertex* vtxESD = fVtxTPC ? esd->GetPrimaryVertexTPC() : esd->GetPrimaryVertexSPD();
  if(!vtxESD)return;
  // Require the vertex to have been reconstructed successfully
  if (strcmp(vtxESD->GetName(), "default")==0)return;
  // vertex position and uncertainties
  fhQA[kVtxPosX] [index]->Fill(vtxESD->GetXv());
  fhQA[kVtxPosY] [index]->Fill(vtxESD->GetYv());
  fhQA[kVtxPosZ] [index]->Fill(vtxESD->GetZv());
  fhQA[kVtxResX] [index]->Fill(vtxESD->GetXRes());
  fhQA[kVtxResY] [index]->Fill(vtxESD->GetYRes());
  fhQA[kVtxResZ] [index]->Fill(vtxESD->GetZRes());
  fhQA[kVtxNCtrb][index]->Fill(vtxESD->GetNContributors());
  
}

//____________________________________________________________________
void AliCFEventRecCuts::SetHistogramBins(Int_t index, Int_t nbins, Double_t *bins)
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
void AliCFEventRecCuts::SetHistogramBins(Int_t index, Int_t nbins, Double_t xmin, Double_t xmax)
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
 void AliCFEventRecCuts::DefineHistograms() {
  //
  // histograms for cut variables
  //
  Int_t color = 2;

  if(!fIsQAOn) {
    AliInfo(Form("No QA histos requested, Please first set the QA flag on!"));
    return;
  }  
  
  // book QA histograms

  Char_t str[256];
  for (Int_t i=0; i<kNStepQA; i++) {
    if (i==0) sprintf(str," ");
    else sprintf(str,"_cut");

    fhQA[kNTracks][i]	= new  TH1F(Form("%s_NTracks%s",GetName(),str),	                "",501,-0.5,500.5);
    fhQA[kVtxPosX][i]	= new  TH1F(Form("%s_Vtx_Pos_X%s",GetName(),str),		"",100,-5.,5.);
    fhQA[kVtxPosY][i]	= new  TH1F(Form("%s_Vtx_Pos_Y%s",GetName(),str),		"",100,-5.,5.);
    fhQA[kVtxPosZ][i]	= new  TH1F(Form("%s_Vtx_Pos_Z%s",GetName(),str),		"",200,-50.,50.);
    fhQA[kVtxResX][i]	= new  TH1F(Form("%s_Vtx_Res_X%s",GetName(),str),		"",100,-1.,1.);
    fhQA[kVtxResY][i]	= new  TH1F(Form("%s_Vtx_Res_Y%s",GetName(),str),		"",100,-1.,1.);
    fhQA[kVtxResZ][i]	= new  TH1F(Form("%s_Vtx_Res_Z%s",GetName(),str),		"",100,-1.,1.);
    fhQA[kVtxNCtrb][i]	= new  TH1F(Form("%s_Vtx_N_Ctrb%s",GetName(),str),		"",1000,0.,1000.);
 
    fhQA[kNTracks][i]	->SetXTitle("Number of ESD tracks");
    fhQA[kVtxPosX][i]	->SetXTitle("Vertex Position X (cm)");
    fhQA[kVtxPosY][i]	->SetXTitle("Vertex Position Y (cm)");
    fhQA[kVtxPosZ][i]	->SetXTitle("Vertex Position Z (cm)");
    fhQA[kVtxResX][i]	->SetXTitle("Vertex Resolution X (cm)");
    fhQA[kVtxResY][i]	->SetXTitle("Vertex Resolution Y (cm)");
    fhQA[kVtxResZ][i]	->SetXTitle("Vertex Resolution Z (cm)");
    fhQA[kVtxNCtrb][i]	->SetXTitle("Number of contributors");
  }

  for(Int_t i=0; i<kNCuts; i++) fhQA[i][1]->SetLineColor(color);

}

//_____________________________________________________________________________
void AliCFEventRecCuts::AddQAHistograms(TList *qaList) {
  //
  // saves the histograms in a TList
  //

  DefineHistograms();

  for (Int_t j=0; j<kNStepQA; j++) {
    for(Int_t i=0; i<kNCuts; i++)
	qaList->Add(fhQA[i][j]);
  }
}
