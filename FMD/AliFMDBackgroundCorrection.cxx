/**************************************************************************
 * Copyright(c) 2004, ALICE Experiment at CERN, All rights reserved. *
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

// Thil class computes background corrections for the FMD. The correction is computed 
// in eta,phi cells and the objects stored can be put into alien to use with analysis.
// It is based on the AliFMDInput class that is used to loop over hits and primaries.
//
// Author: Hans Hjersing Dalsgaard, NBI, hans.dalsgaard@cern.ch
//
//

#include "AliSimulation.h"
#include "TStopwatch.h"
#include "iostream"
#include "TGrid.h"
#include "AliRunLoader.h"
#include "AliGeomManager.h"
#include "AliFMDGeometry.h"
#include "AliStack.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "AliRun.h"
#include "AliFMDBackgroundCorrection.h"
#include "TSystem.h"
#include "AliCDBManager.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TBranch.h"
#include "AliFMDHit.h"
#include "AliLoader.h" 
#include "AliFMD.h"
#include "AliGenEventHeader.h"
#include "AliHeader.h"
#include "TFile.h"
#include "TAxis.h"
#include "AliCDBId.h"
#include "AliCDBMetaData.h"
#include "TROOT.h"
#include "AliFMDParameters.h"
#include "AliLog.h"
#include "TList.h"
#include "AliFMDAnaParameters.h"
#include "AliFMDAnaCalibBackgroundCorrection.h"
#include "AliTrackReference.h"
#include "AliFMDStripIndex.h"

ClassImp(AliFMDBackgroundCorrection)
//_____________________________________________________________________
AliFMDBackgroundCorrection::AliFMDBackgroundCorrection() : 
  TNamed(),
  fCorrectionArray(),
  fPrimaryList()
{} 

//_____________________________________________________________________
AliFMDBackgroundCorrection::AliFMDInputBG::AliFMDInputBG(Bool_t hits_not_trackref) : 
  AliFMDInput("galice.root"),
  fPrimaryArray(),
  fHitArray(),
  fHitMap(),
  fPrim(0),
  fHits(0),
  fZvtxCut(0),
  fNvtxBins(0),
  fPrevTrack(-1),
  fPrevDetector(-1),
  fPrevRing('Q'),
  fPrevSec(-1),
  fNbinsEta(100)
{
  AddLoad(kTracks); 
  if(hits_not_trackref)
    AddLoad(kHits);
  else
    AddLoad(kTrackRefs);
  AddLoad(kKinematics); 
  AddLoad(kHeader);
  
  
}

//_____________________________________________________________________

void AliFMDBackgroundCorrection::GenerateBackgroundCorrection(Bool_t from_hits,
							      const Int_t nvtxbins,
							      Float_t zvtxcut, 
							      const Int_t nBinsEta, 
							      Bool_t storeInAlien, 
							      Int_t runNo,
							      Int_t endRunNo, 
							      const Char_t* filename, 
							      Bool_t simulate,
							      Int_t nEvents,
							      Bool_t inFile,
							      const Char_t* infilename) {
  
  //TGrid::Connect("alien:",0,0,"t");
  if(simulate)
    Simulate(nEvents);
  else {
    //AliCDBManager::Instance()->SetDefaultStorage("alien://Folder=/alice/data/2008/LHC08d/OCDB/");
    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
    AliCDBManager::Instance()->SetRun(runNo);
    
#if defined(__CINT__)
    gSystem->Load("liblhapdf");
    gSystem->Load("libEGPythia6");
    gSystem->Load("libpythia6");
    gSystem->Load("libAliPythia6");
    gSystem->Load("libgeant321");
#endif
      // 
  }  
  
  //Setting up the geometry
  //-----------------------------------------------
  if (AliGeomManager::GetGeometry() == NULL)
    AliGeomManager::LoadGeometry();
  
  AliFMDGeometry* geo = AliFMDGeometry::Instance();
  geo->Init();
  geo->InitTransformations();
    
  AliInfo("Processing hits and primaries ");
  
  AliFMDInputBG input(from_hits);
  
  if(!inFile) {
  
    input.SetVtxCutZ(zvtxcut);
    input.SetNvtxBins(nvtxbins);
    input.SetNbinsEta(nBinsEta);
    input.Run();
  }
  
  AliInfo(Form("Found %d primaries and %d hits.", input.GetNprim(),input.GetNhits()));
  TObjArray* hitArray ;
  TObjArray* primaryArray;
  if(inFile) {
    TFile* infile = TFile::Open(infilename);
    hitArray     = new TObjArray();
    primaryArray = new TObjArray();
    
    for(Int_t det =1; det<=3;det++)
      {
	TObjArray* detArrayHits = new TObjArray();
	detArrayHits->SetName(Form("FMD%d",det));
	hitArray->AddAtAndExpand(detArrayHits,det);
	Int_t nRings = (det==1 ? 1 : 2);
	for(Int_t ring = 0;ring<nRings;ring++)
	  {
	    
	    Char_t ringChar = (ring == 0 ? 'I' : 'O');
	    TObjArray* vtxArrayHits = new TObjArray();
	    vtxArrayHits->SetName(Form("FMD%d%c",det,ringChar));
	    detArrayHits->AddAtAndExpand(vtxArrayHits,ring);
	    for(Int_t v=0; v<nvtxbins;v++)
	      {
	      
		TH2F* hHits          = (TH2F*)infile->Get(Form("hHits_FMD%d%c_vtx%d",det,ringChar,v));
		
		
	      vtxArrayHits->AddAtAndExpand(hHits,v);
	      	      
	    } 
	}
      
    }
    
    for(Int_t iring = 0; iring<2;iring++) {
      Char_t ringChar = (iring == 0 ? 'I' : 'O');
      TObjArray* ringArray = new TObjArray();
      ringArray->SetName(Form("FMD_%c",ringChar));
      primaryArray->AddAtAndExpand(ringArray,iring);
      for(Int_t v=0; v<nvtxbins;v++) {
	
	TH2F* hPrimary       = (TH2F*)infile->Get(Form("hPrimary_FMD_%c_vtx%d",ringChar,v));
	ringArray->AddAtAndExpand(hPrimary,v);
      }
    }
    
    
  }
  else {
    hitArray     = input.GetHits();
    primaryArray = input.GetPrimaries();
  }
  fCorrectionArray.SetName("FMD_bg_correction");
  fCorrectionArray.SetOwner();
   
  TList* primaryList     = new TList();
  primaryList->SetName("primaries");
  
  TList* hitList     = new TList();
  hitList->SetName("hits");
  TList* corrList    = new TList();
  corrList->SetName("corrections");
  
  AliFMDAnaCalibBackgroundCorrection* background = new AliFMDAnaCalibBackgroundCorrection();
  
  for(Int_t det= 1; det <=3; det++) {
    Int_t nRings = (det==1 ? 1 : 2);
    
    TObjArray* detArrayCorrection = new TObjArray();
    detArrayCorrection->SetName(Form("FMD%d",det));
    fCorrectionArray.AddAtAndExpand(detArrayCorrection,det);
    
    
    for(Int_t iring = 0; iring<nRings; iring++) {
      TObjArray* primRingArray = (TObjArray*)primaryArray->At(iring);
      Char_t ringChar = (iring == 0 ? 'I' : 'O');
      TObjArray* vtxArrayCorrection = new TObjArray();
      vtxArrayCorrection->SetName(Form("FMD%d%c",det,ringChar));
      detArrayCorrection->AddAtAndExpand(vtxArrayCorrection,iring);
      
      for(Int_t vertexBin=0;vertexBin<nvtxbins;vertexBin++) {
	TObjArray* detArray  = (TObjArray*)hitArray->At(det);
	TObjArray* vtxArray  = (TObjArray*)detArray->At(iring);
	TH2F* hHits          = (TH2F*)vtxArray->At(vertexBin);
	hitList->Add(hHits);
	TH2F* hPrimary  = (TH2F*)primRingArray->At(vertexBin);
	primaryList->Add(hPrimary);
	TH2F* hCorrection = (TH2F*)hHits->Clone(Form("FMD%d%c_vtxbin_%d_correction",det,ringChar,vertexBin));
	hCorrection->Divide(hPrimary);
	corrList->Add(hCorrection);
	hCorrection->SetTitle(hCorrection->GetName());
	vtxArrayCorrection->AddAtAndExpand(hCorrection,vertexBin);
	background->SetBgCorrection(det,ringChar,vertexBin,hCorrection);
      }
      
    }
  }
  
  TAxis refAxis(nvtxbins,-1*zvtxcut,zvtxcut);
  if(inFile) {
    TFile* infile = TFile::Open(infilename);
    TAxis* refaxis = (TAxis*)infile->Get("vertexbins");
    background->SetRefAxis(refaxis);
      
  }
  else background->SetRefAxis(&refAxis);
  
  TFile*  fout = new TFile(filename,"RECREATE");
  refAxis.Write("vertexbins");
  
  hitList->Write();
  primaryList->Write();
  corrList->Write();
   
  TObjArray* container = new TObjArray();
  container->SetOwner();
  container->AddAtAndExpand(&refAxis,0);
  container->AddAtAndExpand(&fCorrectionArray,1);
  container->AddAtAndExpand(hitArray,2);
  container->AddAtAndExpand(primaryArray,3);
  
  
  if(storeInAlien) {
    AliCDBManager* cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local://$ALICE_ROOT");
    AliCDBId      id(AliFMDAnaParameters::GetBackgroundPath(),runNo,endRunNo);
    
    AliCDBMetaData* meta = new AliCDBMetaData;				
    meta->SetResponsible(gSystem->GetUserInfo()->fRealName.Data());	
    meta->SetAliRootVersion(gROOT->GetVersion());			
    meta->SetBeamPeriod(1);						
    meta->SetComment("Background Correction for FMD");
    meta->SetProperty("key1", background );
    cdb->Put(background, id, meta);
    
  }
  
  fout->Close();
 
  
  }

//_____________________________________________________________________
void AliFMDBackgroundCorrection::Simulate(Int_t nEvents) {
  
  AliSimulation sim ; 
  sim.SetRunNumber(0);
  TGrid::Connect("alien:",0,0,"t");
  sim.SetDefaultStorage("alien://Folder=/alice/data/2008/LHC08d/OCDB/");
  sim.SetConfigFile("Config.C");
  sim.SetRunQA("FMD:");
  TStopwatch timer;
  timer.Start();
  sim.RunSimulation(nEvents);    
  timer.Stop();
  timer.Print();
  
}

//_____________________________________________________________________
Bool_t 
AliFMDBackgroundCorrection::AliFMDInputBG::ProcessHit(AliFMDHit* h, 
						      TParticle* /*p*/) 
{
  
  if(!h)
    return kTRUE;
  Bool_t retval = ProcessEvent(h->Detector(),
			       h->Ring(),
			       h->Sector(),
			       h->Strip(),
			       h->Track(),
			       h->Q());
  
  return retval;
}
//_____________________________________________________________________
Bool_t 
AliFMDBackgroundCorrection::AliFMDInputBG::ProcessTrackRef(AliTrackReference* tr, 
							   TParticle* p) 
{
  if(!tr)
    return kTRUE;
  UShort_t det,sec,strip;
  Char_t   ring;
  AliFMDStripIndex::Unpack(tr->UserId(),det,ring,sec,strip);
  Int_t         nTrack  = tr->GetTrack();
  TDatabasePDG* pdgDB   = TDatabasePDG::Instance();
  TParticlePDG* pdgPart = pdgDB->GetParticle(p->GetPdgCode());
  Float_t       charge  = (pdgPart ? pdgPart->Charge() : 0);
  Bool_t        retval  = ProcessEvent(det,ring,sec,strip,nTrack,charge);
  return retval;
  
}
//_____________________________________________________________________
Bool_t 
AliFMDBackgroundCorrection::AliFMDInputBG::ProcessEvent(UShort_t det,
							Char_t   ring, 
							UShort_t sec, 
							UShort_t strip,
							Int_t    nTrack,
							Float_t  charge)
{
  
  
  
  if(charge !=  0 && 
     ((nTrack != fPrevTrack) || 
      (det != fPrevDetector) || 
      (ring != fPrevRing)    ||
      (sec != fPrevSec))) {
    fHitMap.operator()(det,ring,sec,strip) = 1;
    fHits++;
    
  }
  
  fPrevDetector = det;
  fPrevRing     = ring;
  fPrevSec      = sec;
  fPrevTrack    = nTrack;
  
  return kTRUE;

}

//_____________________________________________________________________
Bool_t AliFMDBackgroundCorrection::AliFMDInputBG::Init() 
{
  fPrimaryArray.SetOwner();
  fPrimaryArray.SetName("FMD_primary");
  
  fPrimaryMapInner.SetBins(fNbinsEta, -6,6, 20, 0, 2*TMath::Pi());
  fPrimaryMapOuter.SetBins(fNbinsEta, -6,6, 40, 0, 2*TMath::Pi());
  fPrimaryMapInner.SetName("fPrimaryMapInner");
  fPrimaryMapInner.SetName("fPrimaryMapOuter");
  
  fPrimaryMapInner.Sumw2();
  fPrimaryMapOuter.Sumw2();
  for(Int_t iring = 0; iring<2;iring++) {
    Char_t ringChar = (iring == 0 ? 'I' : 'O');
    TObjArray* ringArray = new TObjArray();
    ringArray->SetName(Form("FMD_%c",ringChar));
    fPrimaryArray.AddAtAndExpand(ringArray,iring);
    Int_t nSec = (iring == 1 ? 40 : 20);
    for(Int_t v=0; v<fNvtxBins;v++) {

      TH2F* hPrimary       = new TH2F(Form("hPrimary_FMD_%c_vtx%d",ringChar,v),
				      Form("hPrimary_FMD_%c_vtx%d",ringChar,v),
				      fNbinsEta, -6,6, nSec, 0,2*TMath::Pi());
      hPrimary->Sumw2();
      ringArray->AddAtAndExpand(hPrimary,v);
    }
  }
  
  
  fHitArray.SetOwner();
  fHitArray.SetName("FMD_hits");
   
  for(Int_t det =1; det<=3;det++) {
    TObjArray* detArrayHits = new TObjArray();
    detArrayHits->SetName(Form("FMD%d",det));
    fHitArray.AddAtAndExpand(detArrayHits,det);
    Int_t nRings = (det==1 ? 1 : 2);
    for(Int_t ring = 0;ring<nRings;ring++) {
      Int_t nSec = (ring == 1 ? 40 : 20);
      Char_t ringChar = (ring == 0 ? 'I' : 'O');
      TObjArray* vtxArrayHits = new TObjArray();
      vtxArrayHits->SetName(Form("FMD%d%c",det,ringChar));
      detArrayHits->AddAtAndExpand(vtxArrayHits,ring);
      for(Int_t v=0; v<fNvtxBins;v++) {
	TH2F* hHits = new TH2F(Form("hHits_FMD%d%c_vtx%d", det,ringChar,v),
			       Form("hHits_FMD%d%c_vtx%d", det,ringChar,v),
			       fNbinsEta, -6,6, nSec, 0, 2*TMath::Pi());
	hHits->Sumw2();
	vtxArrayHits->AddAtAndExpand(hHits,v);
	
      } 
    }
  }

  AliFMDInput::Init();
  
  return kTRUE;
}
//_____________________________________________________________________

Bool_t AliFMDBackgroundCorrection::AliFMDInputBG::Begin(Int_t event ) 
{

  Bool_t             retVal    = AliFMDInput::Begin(event); 
  AliStack*          partStack = fLoader->Stack();
  Int_t              nTracks   = partStack->GetNtrack();
  
  
  for(Int_t j=0;j<nTracks;j++) {
    TParticle* p           = partStack->Particle(j);
    TDatabasePDG* pdgDB    = TDatabasePDG::Instance();
    TParticlePDG* pdgPart  = pdgDB->GetParticle(p->GetPdgCode());
    Float_t       charge   = (pdgPart ? pdgPart->Charge() : 0);
    Float_t       phi      = TMath::ATan2(p->Py(),p->Px());
    
    if(phi<0) phi = phi+2*TMath::Pi();
    // if(p->Theta() == 0) continue;
    Float_t eta   = p->Eta();   
    
    // std::cout<<-1*TMath::Log(TMath::Tan(0.5*p->Theta()))<<std::endl;
    
    Bool_t primary = partStack->IsPhysicalPrimary(j);
    //(charge!=0)&&(TMath::Abs(p->Vx() - vertex.At(0))<0.01)&&(TMath::Abs(p->Vy() - vertex.At(1))<0.01)&&(TMath::Abs(p->Vz() - vertex.At(2))<0.01);
    if(primary && charge !=0) {
      
      
      
      fPrim++;
      
      fPrimaryMapInner.Fill(eta,phi);
      fPrimaryMapOuter.Fill(eta,phi);
    }
  }
  
  return retVal;
}
//_____________________________________________________________________
Bool_t AliFMDBackgroundCorrection::AliFMDInputBG::End()  {
  
  Bool_t retval = AliFMDInput::End();
  
  AliGenEventHeader* genHeader = fLoader->GetHeader()->GenEventHeader();
  TArrayF vertex;
  genHeader->PrimaryVertex(vertex);
  
  if(TMath::Abs(vertex.At(2)) > fZvtxCut) 
    return kTRUE;
  
  Double_t delta           = 2*fZvtxCut/fNvtxBins;
  Double_t vertexBinDouble = (vertex.At(2) + fZvtxCut) / delta;
  Int_t    vertexBin       = (Int_t)vertexBinDouble;
  //Primaries
  TObjArray* innerArray = (TObjArray*)fPrimaryArray.At(0);
  TObjArray* outerArray = (TObjArray*)fPrimaryArray.At(1);
  
  TH2F* hPrimaryInner  = (TH2F*)innerArray->At(vertexBin);
  TH2F* hPrimaryOuter  = (TH2F*)outerArray->At(vertexBin);
  
  hPrimaryInner->Add(&fPrimaryMapInner);
  hPrimaryOuter->Add(&fPrimaryMapOuter);
  
  //Hits
  for(UShort_t det=1;det<=3;det++) {
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t   ring = (ir == 0 ? 'I' : 'O');
      UShort_t nsec = (ir == 0 ? 20  : 40);
      UShort_t nstr = (ir == 0 ? 512 : 256);
      
      for(UShort_t sec =0; sec < nsec;  sec++) {
	
	for(UShort_t strip = 0; strip < nstr; strip++) {
	  
	  if(fHitMap.operator()(det,ring,sec,strip) > 0) {
	  
	    Double_t x,y,z;
	    AliFMDGeometry* fmdgeo = AliFMDGeometry::Instance();
	    fmdgeo->Detector2XYZ(det,ring,sec,strip,x,y,z);
	    
	    Int_t iring = (ring == 'I' ? 0 : 1);
	    
	    TObjArray* detArray  = (TObjArray*)fHitArray.At(det);
	    TObjArray* vtxArray  = (TObjArray*)detArray->At(iring);
	    TH2F* hHits          = (TH2F*)vtxArray->At(vertexBin);
	    
	    Float_t   phi   = TMath::ATan2(y,x);
	    if(phi<0) phi   = phi+2*TMath::Pi();
	    Float_t   r     = TMath::Sqrt(TMath::Power(x,2)+TMath::Power(y,2));
	    Float_t   theta = TMath::ATan2(r,z-vertex.At(2));
	    Float_t   eta   = -1*TMath::Log(TMath::Tan(0.5*theta));
	    hHits->Fill(eta,phi);
	  }
	  
	}
      }
    }
  }
  
  fPrimaryMapInner.Reset();
  fPrimaryMapOuter.Reset();
  fHitMap.Reset(0);
  
  return retval;
}

//
// EOF
//


