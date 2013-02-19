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
#include "AliAODEvent.h"
#include "AliAODVertex.h"
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
  fVtxSPD(0),
  fCheckCorrelationSPDVtx(0),
  fVtxResolution(kFALSE),
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
  fVtxSPD(0),
  fCheckCorrelationSPDVtx(0),
  fVtxResolution(kFALSE),
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
  fVtxSPD(c.fVtxSPD),
  fCheckCorrelationSPDVtx(c.fCheckCorrelationSPDVtx),
  fVtxResolution(c.fVtxResolution),
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
    fVtxSPD=c.fVtxSPD;
    fCheckCorrelationSPDVtx=c.fCheckCorrelationSPDVtx;
    fVtxResolution = c.fVtxResolution;
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

  AliVEvent *inputEvent = dynamic_cast<AliVEvent *>(obj);
  if(!inputEvent){
    AliDebug(1, "Not a virtual event");
    for(Int_t j=0;j<kNCuts;j++)fBitMap->SetBitNumber(j,kFALSE);
    return;
  }
  const AliVVertex *vtxTracks = GetPrimaryVertexTracks(inputEvent),
                   *vtxSPD = GetPrimaryVertexSPD(inputEvent),
                   *vtxPrim(NULL);
  
  //first assume the event will be accepted: 
  for(Int_t j=0;j<kNCuts;j++)fBitMap->SetBitNumber(j,kTRUE);

  if(fVtxMixed){
    // Use mixed vertex: Prefer vertex with tracks, in case not available use SPD vertex
    if(vtxTracks && vtxTracks->GetNContributors() > 0) vtxPrim = vtxTracks;
    else if(vtxSPD && vtxSPD->GetNContributors() > 0) vtxPrim = vtxSPD;
  } else if(fVtxSPD){
    if(vtxSPD && vtxSPD->GetNContributors () > 0) vtxPrim = vtxSPD;
  } else {
    if(vtxTracks && vtxTracks->GetNContributors() > 0) vtxPrim = vtxTracks;
  }
  if(!vtxPrim){
    // No primary vertex: Reject event
	  for(Int_t j=1;j<kNCuts;j++)fBitMap->SetBitNumber(j,kFALSE); 
	  AliWarning("Cannot get vertex, skipping event");
	  return;
  }
  // Check standard vertex cuts using the primary vertex 
  if(fVtxZMin > - 1000. && fVtxZMax < 1000.){
    // Primary vertex z cut required
    // Pick up the position and uncertainties
    Double_t vtxPos[3];
    vtxPos[0] = vtxPrim->GetX();
    vtxPos[1] = vtxPrim->GetY();
    vtxPos[2] = vtxPrim->GetZ();
    if (vtxPos[2]>fVtxZMax || vtxPos[2]<fVtxZMin)
	    fBitMap->SetBitNumber(kVtxPosZ,kFALSE); 
  }
  if(fVtxNCtrbMin > 0){  
    // cut required if the cut value is set to something larger than 0
    // same effect as setting the min. cut value to 0
    Int_t nCtrb = vtxPrim->GetNContributors();
    if (nCtrb<fVtxNCtrbMin)
	    fBitMap->SetBitNumber(kVtxNCtrb,kFALSE);
  }

  // check vertex correlation cut
  if(fCheckCorrelationSPDVtx){
    if(vtxTracks && vtxTracks->GetNContributors() && vtxSPD && vtxSPD->GetNContributors()){
      // Both vertices available, check correlation
      if(TMath::Abs(vtxTracks->GetZ() - vtxSPD->GetZ()) >= 0.5) fBitMap->SetBitNumber(kCorrelation, kFALSE);
    } else {
      // No correlation available: set the cut to false
      fBitMap->SetBitNumber(kCorrelation, kFALSE); 
    }
  }
  // check cut on the SPD vertex resolution
  if(fVtxResolution){
    if(vtxSPD){
      TString vtxTyp = vtxSPD->GetTitle();
      Double_t cov[6]={0};
      vtxSPD->GetCovarianceMatrix(cov);
      Double_t zRes = TMath::Sqrt(cov[5]);
      if (vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) fBitMap->SetBitNumber(kResolution, kFALSE);
    } else{
      fBitMap->SetBitNumber(kResolution, kFALSE);
    }
  }
}

//_____________________________________________________________________________
void AliHFEextraEventCuts::FillHistograms(TObject* obj, Bool_t b)
{
  //
  // fill the QA histograms
  //

  if(!fIsQAOn) return;
  // cast TObject into VEvent
  AliVEvent* inputEvent = dynamic_cast<AliESDEvent *>(obj);
  if (!inputEvent) return;

  // index = 0: fill histograms before cuts
  // index = 1: fill histograms after cuts
  Int_t index = -1;
  index = ((b) ? 1 : 0);

  // Obtain vertices
  const AliVVertex *vtxTracks = GetPrimaryVertexTracks(inputEvent),
                   *vtxSPD = GetPrimaryVertexSPD(inputEvent),
                   *vtxPrim(NULL);

  //look at vertex parameters:
  if(fVtxMixed){
    if(vtxTracks && vtxTracks->GetNContributors() > 0) vtxPrim = vtxTracks;
    else if(vtxSPD && vtxSPD->GetNContributors()) vtxPrim = vtxSPD;
  }
  else {   
    vtxPrim = vtxTracks; 
  }
  if(!vtxPrim)return;
  // vertex position and uncertainties
  fhQA[kVtxPosZ] [index]->Fill(vtxPrim->GetZ());
  fhQA[kVtxNCtrb][index]->Fill(vtxPrim->GetNContributors());
  // SPD Vertex Position Correlation 
  if(vtxTracks && vtxSPD){
    fhQA[kCorrelation][index]->Fill(vtxTracks->GetZ()-vtxSPD->GetZ());
  }
  if(vtxSPD){
    Double_t cov[6]={0};
    vtxSPD->GetCovarianceMatrix(cov);
    fhQA[kResolution][index]->Fill(TMath::Sqrt(cov[5]));
  }
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
    fhQA[kCorrelation][i] = new TH1F(Form("%s_SPDCorrelation_%s",GetName(),str), "",200, -10., 10);
    fhQA[kResolution][i] = new TH1F(Form("%s_SPDResolution_%s",GetName(), str), "", 100, 0., 1.); 
 
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

//_____________________________________________________________________________
const AliVVertex *AliHFEextraEventCuts::GetPrimaryVertexSPD(const AliVEvent * const inputEvent){
  //
  // Get SPD vertex from event
  //
  const AliVVertex *spdvtx(NULL);
  const AliESDEvent *esd(NULL);
  const AliAODEvent *aod(NULL);
  if((esd = dynamic_cast<const AliESDEvent *>(inputEvent))){
    spdvtx = esd->GetPrimaryVertexSPD();
  } else if((aod = dynamic_cast<const AliAODEvent *>(inputEvent))){
    spdvtx = aod->GetPrimaryVertexSPD();
  }
  return spdvtx;
}

//_____________________________________________________________________________
const AliVVertex *AliHFEextraEventCuts::GetPrimaryVertexTracks(const AliVEvent *const inputEvent){
  // 
  // Get Primary Vertex from tracks
  //
  const AliVVertex *trkvtx(NULL);
  const AliESDEvent *esd(NULL);
  const AliAODEvent *aod(NULL);
  if((esd = dynamic_cast<const AliESDEvent *>(inputEvent))){
    trkvtx = esd->GetPrimaryVertexTracks();
  } else if((aod = dynamic_cast<const AliAODEvent *>(inputEvent))){
    const AliVVertex *vtxTmp = aod->GetPrimaryVertex();
    // check whether the primary vertex is the vertex from tracks
    TString vtxTtl = vtxTmp->GetTitle();
    if(vtxTtl.Contains("VertexerTracks")){
      trkvtx = vtxTmp;
    }
  }
  return trkvtx;
}

