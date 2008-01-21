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
  fBitMap(0x0),
  fhNBinsNTracks(0),
  fhBinLimNTracks(0),
  fhNBinsVtxPosX(0),
  fhBinLimVtxPosX(0),
  fhNBinsVtxPosY(0),
  fhBinLimVtxPosY(0),
  fhNBinsVtxPosZ(0),
  fhBinLimVtxPosZ(0),
  fhNBinsVtxResX(0),
  fhBinLimVtxResX(0),
  fhNBinsVtxResY(0),
  fhBinLimVtxResY(0),
  fhNBinsVtxResZ(0),
  fhBinLimVtxResZ(0)
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
  fBitMap(0x0),
  fhNBinsNTracks(0),
  fhBinLimNTracks(0),
  fhNBinsVtxPosX(0),
  fhBinLimVtxPosX(0),
  fhNBinsVtxPosY(0),
  fhBinLimVtxPosY(0),
  fhNBinsVtxPosZ(0),
  fhBinLimVtxPosZ(0),
  fhNBinsVtxResX(0),
  fhBinLimVtxResX(0),
  fhNBinsVtxResY(0),
  fhBinLimVtxResY(0),
  fhNBinsVtxResZ(0),
  fhBinLimVtxResZ(0)
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
  fBitMap(c.fBitMap),
  fhNBinsNTracks(c.fhNBinsNTracks),
  fhBinLimNTracks(c.fhBinLimNTracks),
  fhNBinsVtxPosX(c.fhNBinsVtxPosX),
  fhBinLimVtxPosX(c.fhBinLimVtxPosX),
  fhNBinsVtxPosY(c.fhNBinsVtxPosY),
  fhBinLimVtxPosY(c.fhBinLimVtxPosY),
  fhNBinsVtxPosZ(c.fhNBinsVtxPosZ),
  fhBinLimVtxPosZ(c.fhBinLimVtxPosZ),
  fhNBinsVtxResX(c.fhNBinsVtxResX),
  fhBinLimVtxResX(c.fhBinLimVtxResX),
  fhNBinsVtxResY(c.fhNBinsVtxResY),
  fhBinLimVtxResY(c.fhBinLimVtxResY),
  fhNBinsVtxResZ(c.fhNBinsVtxResZ),
  fhBinLimVtxResZ(c.fhBinLimVtxResZ) 
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

  if(fhBinLimNTracks)delete fhBinLimNTracks;
  if(fhBinLimVtxPosX)delete fhBinLimVtxPosX;
  if(fhBinLimVtxPosY)delete fhBinLimVtxPosY;
  if(fhBinLimVtxPosZ)delete fhBinLimVtxPosZ;
  if(fhBinLimVtxResX)delete fhBinLimVtxResX;
  if(fhBinLimVtxResY)delete fhBinLimVtxResY;
  if(fhBinLimVtxResZ)delete fhBinLimVtxResZ;

  if(fBitMap)delete fBitMap;

}

//_____________________________________________________________________________
void AliCFEventRecCuts::Init() {
  //
  // initialises all QA histograms
  //
  if(fIsQAOn)
    DefineHistograms();
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

  //set default bin number/ranges for QA histograms

  SetHistogramBins(kNTracks,23,-0.5,22.5);
  SetHistogramBins(kVtxPosX,100,-5,5);
  SetHistogramBins(kVtxPosY,100,-5,5);
  SetHistogramBins(kVtxPosZ,100,-50,50);
  SetHistogramBins(kVtxResX,100,-1,1);
  SetHistogramBins(kVtxResY,100,-1,1);
  SetHistogramBins(kVtxResZ,100,-1,1);

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
    fBitMap=c.fBitMap;
    fhNBinsNTracks=c.fhNBinsNTracks;
    fhBinLimNTracks=c.fhBinLimNTracks;
    fhNBinsVtxPosX=c.fhNBinsVtxPosX;
    fhBinLimVtxPosX=c.fhBinLimVtxPosX;
    fhNBinsVtxPosY=c.fhNBinsVtxPosY;
    fhBinLimVtxPosY=c.fhBinLimVtxPosY;
    fhNBinsVtxPosZ=c.fhNBinsVtxPosZ;
    fhBinLimVtxPosZ=c.fhBinLimVtxPosZ;
    fhNBinsVtxResX=c.fhNBinsVtxResX;
    fhBinLimVtxResX=c.fhBinLimVtxResX;
    fhNBinsVtxResY=c.fhNBinsVtxResY;
    fhBinLimVtxResY=c.fhBinLimVtxResY;
    fhNBinsVtxResZ=c.fhNBinsVtxResZ;
    fhBinLimVtxResZ=c.fhBinLimVtxResZ;
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

  TBits *bitmap = SelectionBitMap(obj);

  Bool_t isSelected = kTRUE;

  for (UInt_t icut=0; icut<bitmap->GetNbits();icut++)
	if(!bitmap->TestBitNumber(icut)) isSelected = kFALSE;

  return isSelected;

}

//____________________________________________________________________
TBits *AliCFEventRecCuts::SelectionBitMap(TObject* obj) {
  //
  //cut on the number of charged tracks and on the event vertex.
  //so far specific to AliESDEvents
  //

  //Check if the requested cuts are passed and return a bitmap
  for(Int_t j=0;j<kNCuts;j++)fBitMap->SetBitNumber(j,kFALSE);
  AliESDEvent* esd = dynamic_cast<AliESDEvent *>(obj);
  if ( !esd ) return fBitMap ;

  //now start checking the cuts,
  //first assume the event will be accepted: 
  for(Int_t j=0;j<kNCuts;j++)fBitMap->SetBitNumber(j,kTRUE);

  //Number of charged tracks:
  Int_t nTracks = esd->GetNumberOfTracks();
  if(nTracks<fNTracksMin || nTracks>fNTracksMax)
    fBitMap->SetBitNumber(0,kFALSE); 
  
  if(fRequireVtxCuts){
    const AliESDVertex* vtxESD = esd->GetVertex();
    if(!vtxESD){
      for(Int_t j=1;j<kNCuts;j++)fBitMap->SetBitNumber(j,kFALSE); 
      return fBitMap;
    }
    // Require the vertex to have been reconstructed successfully
    if (strcmp(vtxESD->GetName(), "default")==0){
      AliWarning(Form(" No reconstructed vertex found, skip event"));    
      for(Int_t j=1;j<kNCuts;j++)fBitMap->SetBitNumber(j,kFALSE); 
      return fBitMap;
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
  }  
  return fBitMap;
}

//_____________________________________________________________________________
void AliCFEventRecCuts::GetBitMap(TObject* obj, TBits *bitmap) {
  //
  // retrieve the pointer to the bitmap
  //
  bitmap = SelectionBitMap(obj);

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
  const AliESDVertex* vtxESD = esd->GetVertex();
  if(!vtxESD)return;
  // Require the vertex to have been reconstructed successfully
  if (strcmp(vtxESD->GetName(), "default")==0)return;
  // vertex position and uncertainties
  fhQA[kVtxPosX][index]->Fill(vtxESD->GetXv());
  fhQA[kVtxPosY][index]->Fill(vtxESD->GetYv());
  fhQA[kVtxPosZ][index]->Fill(vtxESD->GetZv());
  fhQA[kVtxResX][index]->Fill(vtxESD->GetXRes());
  fhQA[kVtxResY][index]->Fill(vtxESD->GetYRes());
  fhQA[kVtxResZ][index]->Fill(vtxESD->GetZRes());
  
}

//_____________________________________________________________________________
void AliCFEventRecCuts::SetHistogramBins(Int_t index, Int_t nbins, Double_t *bins)
{
  //
  // QA histogram axis parameters
  // variable bin size:user inputs nbins and the vector of bin limits
  //

  switch(index){
  case kNTracks:
    fhNBinsNTracks=nbins+1;
    fhBinLimNTracks=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimNTracks[i]=bins[i];
    break;
  case kVtxPosX:
    fhNBinsVtxPosX=nbins+1;
    fhBinLimVtxPosX=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimVtxPosX[i]=bins[i];
    break;
  case kVtxPosY:
    fhNBinsVtxPosY=nbins+1;
    fhBinLimVtxPosY=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimVtxPosY[i]=bins[i];
    break;
  case kVtxPosZ:
    fhNBinsVtxPosZ=nbins+1;
    fhBinLimVtxPosZ=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimVtxPosZ[i]=bins[i];
    break;
  case kVtxResX:
    fhNBinsVtxResX=nbins+1;
    fhBinLimVtxResX=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimVtxResX[i]=bins[i];
    break;
  case kVtxResY:
    fhNBinsVtxResY=nbins+1;
    fhBinLimVtxResY=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimVtxResY[i]=bins[i];
    break;
  case kVtxResZ:
    fhNBinsVtxResZ=nbins+1;
    fhBinLimVtxResZ=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimVtxResZ[i]=bins[i];
    break;
  }

}

//_____________________________________________________________________________
void AliCFEventRecCuts::SetHistogramBins(Int_t index, Int_t nbins, Double_t xmin, Double_t xmax)
{
  //
  // QA histogram axis parameters
  // fixed bin size: user inputs nbins, xmin and xmax
  //
  switch(index){
  case kNTracks:
    fhNBinsNTracks=nbins+1;
    fhBinLimNTracks=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimNTracks[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;
  case kVtxPosX:
    fhNBinsVtxPosX=nbins+1;
    fhBinLimVtxPosX=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimVtxPosX[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;
  case kVtxPosY:
    fhNBinsVtxPosY=nbins+1;
    fhBinLimVtxPosY=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimVtxPosY[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;
  case kVtxPosZ:
    fhNBinsVtxPosZ=nbins+1;
    fhBinLimVtxPosZ=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimVtxPosZ[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;
  case kVtxResX:
    fhNBinsVtxResX=nbins+1;
    fhBinLimVtxResX=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimVtxResX[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;
  case kVtxResY:
    fhNBinsVtxResY=nbins+1;
    fhBinLimVtxResY=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimVtxResY[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;
  case kVtxResZ:
    fhNBinsVtxResZ=nbins+1;
    fhBinLimVtxResZ=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimVtxResZ[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;
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

    fhQA[kNTracks][i]	= new  TH1F(Form("%s_NTracks%s",GetName(),str),	                "",fhNBinsNTracks-1,fhBinLimNTracks);
    fhQA[kVtxPosX][i]	= new  TH1F(Form("%s_Vtx_Pos_X%s",GetName(),str),		"",fhNBinsVtxPosX-1,fhBinLimVtxPosX);
    fhQA[kVtxPosY][i]	= new  TH1F(Form("%s_Vtx_Pos_Y%s",GetName(),str),		"",fhNBinsVtxPosY-1,fhBinLimVtxPosY);
    fhQA[kVtxPosZ][i]	= new  TH1F(Form("%s_Vtx_Pos_Z%s",GetName(),str),		"",fhNBinsVtxPosZ-1,fhBinLimVtxPosZ);

    fhQA[kVtxResX][i]	= new  TH1F(Form("%s_Vtx_Res_X%s",GetName(),str),		"",fhNBinsVtxResX-1,fhBinLimVtxResX);
    fhQA[kVtxResY][i]	= new  TH1F(Form("%s_Vtx_Res_Y%s",GetName(),str),		"",fhNBinsVtxResY-1,fhBinLimVtxResY);
    fhQA[kVtxResZ][i]	= new  TH1F(Form("%s_Vtx_Res_Z%s",GetName(),str),		"",fhNBinsVtxResZ-1,fhBinLimVtxResZ);
 
    fhQA[kNTracks][i]	->SetXTitle("Number of ESD tracks");
    fhQA[kVtxPosX][i]	->SetXTitle("Vertex Position X (cm)");
    fhQA[kVtxPosY][i]	->SetXTitle("Vertex Position Y (cm)");
    fhQA[kVtxPosZ][i]	->SetXTitle("Vertex Position Z (cm)");
    fhQA[kVtxResX][i]	->SetXTitle("Vertex Resolution X (cm)");
    fhQA[kVtxResY][i]	->SetXTitle("Vertex Resolution Y (cm)");
    fhQA[kVtxResZ][i]	->SetXTitle("Vertex Resolution Z (cm)");

  }

  for(Int_t i=0; i<kNCuts; i++) fhQA[i][1]->SetLineColor(color);

}

//_____________________________________________________________________________
void AliCFEventRecCuts::AddQAHistograms(TList *list) const {
  //
  // saves the histograms in a TList
  //
  if(!fIsQAOn) return;  

  for (Int_t j=0; j<kNStepQA; j++) {
    for(Int_t i=0; i<kNCuts; i++)
	list->Add(fhQA[i][j]);
  }
}
