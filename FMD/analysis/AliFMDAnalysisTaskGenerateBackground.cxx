#include "AliFMDAnalysisTaskGenerateBackground.h"
#include "AliESDEvent.h"
#include "iostream"
#include "AliESDFMD.h"
#include "TH2F.h"
#include "AliTrackReference.h"
#include "AliStack.h"
#include "AliFMDAnaParameters.h"
#include "AliFMDStripIndex.h"
#include "AliStack.h"
#include "AliMCParticle.h"
#include "AliMCEvent.h"
//#include "AliFMDGeometry.h"
#include "TArray.h"
#include "AliGenEventHeader.h"
#include "AliHeader.h"
#include "AliFMDAnaCalibBackgroundCorrection.h"
//#include "AliCDBManager.h"
//#include "AliCDBId.h"
//#include "AliCDBMetaData.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TAxis.h"
ClassImp(AliFMDAnalysisTaskGenerateBackground)

//_____________________________________________________________________
AliFMDAnalysisTaskGenerateBackground::AliFMDAnalysisTaskGenerateBackground():
AliAnalysisTaskSE()
{
  // Default constructor
}
//_____________________________________________________________________
AliFMDAnalysisTaskGenerateBackground::AliFMDAnalysisTaskGenerateBackground(const char* name):
  AliAnalysisTaskSE(name),
  fZvtxCut(10),
  fNvtxBins(10),
  fNbinsEta(200)
{
 
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TH1F::Class());
  DefineOutput(4, TList::Class());
}
//_____________________________________________________________________
void AliFMDAnalysisTaskGenerateBackground::UserCreateOutputObjects()
{
// Create the output containers
//
  
  std::cout<<"Creating output objects"<<std::endl;
  for(Int_t iring = 0; iring<2;iring++) {
    Char_t ringChar = (iring == 0 ? 'I' : 'O');
    Int_t nSec = (iring == 1 ? 40 : 20);
    for(Int_t v=0; v<fNvtxBins;v++) {

      TH2F* hPrimary       = new TH2F(Form("hPrimary_FMD_%c_vtx%d",ringChar,v),
				      Form("hPrimary_FMD_%c_vtx%d",ringChar,v),
				      fNbinsEta, -6,6, nSec, 0,2*TMath::Pi());
      hPrimary->Sumw2();
      fListOfPrimaries.Add(hPrimary);
    }
  }
  
  
  for(Int_t det =1; det<=3;det++) {
    Int_t nRings = (det==1 ? 1 : 2);
    for(Int_t ring = 0;ring<nRings;ring++) {
      Int_t nSec = (ring == 1 ? 40 : 20);
      Char_t ringChar = (ring == 0 ? 'I' : 'O');
      for(Int_t v=0; v<fNvtxBins;v++) {
	TH2F* hHits = new TH2F(Form("hHits_FMD%d%c_vtx%d", det,ringChar,v),
			       Form("hHits_FMD%d%c_vtx%d", det,ringChar,v),
			       fNbinsEta, -6,6, nSec, 0, 2*TMath::Pi());
	hHits->Sumw2();
	fListOfHits.Add(hHits);
	
      } 
    }
  }
  TH1F* doubleHits = new TH1F("DoubleHits",  "DoubleHits",
			 fNbinsEta, -6,6);
  TH1F* allHits = new TH1F("allHits",  "allHits",
			 fNbinsEta, -6,6);
  fVertexBins.SetName("VertexBins");
  fVertexBins.GetXaxis()->Set(fNvtxBins,-1*fZvtxCut,fZvtxCut);
  fListOfHits.Add(allHits);
  fListOfHits.Add(doubleHits);
}
//_____________________________________________________________________
void AliFMDAnalysisTaskGenerateBackground::Init()
{
  std::cout<<"Init"<<std::endl;
  
  
  fLastTrackByStrip.Reset(-1);
  
  
}
//_____________________________________________________________________
void AliFMDAnalysisTaskGenerateBackground::UserExec(Option_t */*option*/)
{
  
  fLastTrackByStrip.Reset(-1);
  fHitsByStrip.Reset(0);
  AliMCEvent* mcevent = MCEvent();
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  AliMCParticle* particle = 0;
  AliStack* stack = mcevent->Stack();
  
  UShort_t det,sec,strip;
  Char_t   ring;
  
  Int_t nTracks = mcevent->GetNumberOfTracks();
  AliHeader* header            = mcevent->Header();
  AliGenEventHeader* genHeader = header->GenEventHeader();
  
  TArrayF vertex;
  genHeader->PrimaryVertex(vertex);
  
  if(TMath::Abs(vertex.At(2)) > fZvtxCut)
    return;
  for(Int_t i = 0 ;i<nTracks;i++) {
    particle = mcevent->GetTrack(i);
    
    if(!particle)
      continue;
    
    Double_t delta           = 2*fZvtxCut/fNvtxBins;
    Double_t vertexBinDouble = (vertex.At(2) + fZvtxCut) / delta;
    Int_t    vertexBin       = (Int_t)vertexBinDouble;
    
    if(stack->IsPhysicalPrimary(i) && particle->Charge() != 0) {
      
      
      TH2F* hPrimaryInner = (TH2F*)fListOfPrimaries.FindObject( Form("hPrimary_FMD_%c_vtx%d",'I',vertexBin));
      TH2F* hPrimaryOuter = (TH2F*)fListOfPrimaries.FindObject( Form("hPrimary_FMD_%c_vtx%d",'O',vertexBin));
      hPrimaryInner->Fill(particle->Eta(),particle->Phi());
      hPrimaryOuter->Fill(particle->Eta(),particle->Phi());      
    }
    
    for(Int_t j=0; j<particle->GetNumberOfTrackReferences();j++) {
      
      AliTrackReference* ref = particle->GetTrackReference(j);
      
      if(ref->DetectorId() != AliTrackReference::kFMD)
	continue;
      AliFMDStripIndex::Unpack(ref->UserId(),det,ring,sec,strip);
      Float_t thisStripTrack = fLastTrackByStrip.operator()(det,ring,sec,strip);
      if(particle->Charge() != 0 && i != thisStripTrack ) {
	//Double_t x,y,z;
	//AliFMDGeometry* fmdgeo = AliFMDGeometry::Instance();
	//fmdgeo->Detector2XYZ(det,ring,sec,strip,x,y,z);
	Float_t phi = pars->GetPhiFromSector(det,ring,sec);
	Float_t eta = pars->GetEtaFromStrip(det,ring,sec,strip,vertex.At(2));
	//Float_t   phi   = TMath::ATan2(y,x);
	//if(phi<0) phi   = phi+2*TMath::Pi();
	//	Float_t   r     = TMath::Sqrt(TMath::Power(x,2)+TMath::Power(y,2));
	//Float_t   theta = TMath::ATan2(r,z-vertex.At(2));
	//Float_t   eta   = -1*TMath::Log(TMath::Tan(0.5*theta));
	TH2F* hHits = (TH2F*)fListOfHits.FindObject(Form("hHits_FMD%d%c_vtx%d", det,ring,vertexBin));
	hHits->Fill(eta,phi);
	Float_t nstrips = (ring =='O' ? 256 : 512);
	fHitsByStrip.operator()(det,ring,sec,strip) +=1;
	TH1F* allHits = (TH1F*)fListOfHits.FindObject("allHits");
	TH1F* doubleHits = (TH1F*)fListOfHits.FindObject("DoubleHits");
	
	if(fHitsByStrip.operator()(det,ring,sec,strip) == 1)
	  allHits->Fill(eta);
	
	doubleHits->Fill(eta);
	/*if(fHitsByStrip.operator()(det,ring,sec,strip) == 2){
	  TH1F* doubleHits = (TH1F*)fListOfHits.FindObject("DoubleHits");
	  doubleHits->Fill(eta,2);
	  }*/
	//if(fHitsByStrip.operator()(det,ring,sec,strip) > 1){
	//  doubleHits->Fill(eta);
	//	}
	
	
	fLastTrackByStrip.operator()(det,ring,sec,strip) = (Float_t)i;
	if(strip >0)
	  fLastTrackByStrip.operator()(det,ring,sec,strip-1) = (Float_t)i;
	if(strip < (nstrips - 1))
	  fLastTrackByStrip.operator()(det,ring,sec,strip+1) = (Float_t)i;
      }
    }

  }
    
  	
  PostData(1, &fListOfHits);
  PostData(2, &fListOfPrimaries);
  PostData(3, &fVertexBins);
}
//_____________________________________________________________________
void AliFMDAnalysisTaskGenerateBackground::Terminate(Option_t */*option*/)
{
  /*  TH1F* allHits = (TH1F*)fListOfHits.FindObject("allHits");
  TH1F* doubleHits = (TH1F*)fListOfHits.FindObject("DoubleHits");
  
  doubleHits->Divide(allHits);
  GenerateCorrection();
  PostData(1, &fListOfHits);
  PostData(4, &fListOfCorrection);*/
  
}
//_____________________________________________________________________
void AliFMDAnalysisTaskGenerateBackground::GenerateCorrection() {
  
  fBackground = new AliFMDAnaCalibBackgroundCorrection();
  
  for(Int_t det= 1; det <=3; det++) {
    Int_t nRings = (det==1 ? 1 : 2);
    
    for(Int_t iring = 0; iring<nRings; iring++) {
      Char_t ring = (iring == 0 ? 'I' : 'O');
      for(Int_t vertexBin=0;vertexBin<fNvtxBins  ;vertexBin++) {
	TH2F* hHits          = (TH2F*)fListOfHits.FindObject(Form("hHits_FMD%d%c_vtx%d", det,ring,vertexBin));
	TH2F* hPrimary  = (TH2F*)fListOfPrimaries.FindObject( Form("hPrimary_FMD_%c_vtx%d",ring,vertexBin));
	TH2F* hCorrection = (TH2F*)hHits->Clone(Form("FMD%d%c_vtxbin_%d_correction",det,ring,vertexBin));
	hCorrection->Divide(hPrimary);
	hCorrection->SetTitle(hCorrection->GetName());
	fListOfCorrection.Add(hCorrection);
	fBackground->SetBgCorrection(det,ring,vertexBin,hCorrection);
      }
      
    }
  }
  TAxis refAxis(fNvtxBins,-1*fZvtxCut,fZvtxCut);
  fBackground->SetRefAxis(&refAxis);

}
//_____________________________________________________________________
void AliFMDAnalysisTaskGenerateBackground::ReadFromFile(const Char_t* filename, Bool_t storeInOCDB, Int_t runNo) {

  TFile infile(filename);
  TH1F* hVertex = (TH1F*)infile.Get("VertexBins");
  fZvtxCut = hVertex->GetXaxis()->GetXmax();
  fNvtxBins = hVertex->GetXaxis()->GetNbins();
  fVertexBins.SetName("VertexBins");
  fVertexBins.GetXaxis()->Set(fNvtxBins,-1*fZvtxCut,fZvtxCut);
  
  TList* listOfHits = (TList*)infile.Get("Hits");
  TList* listOfPrim = (TList*)infile.Get("Primaries");
  
  for(Int_t det =1; det<=3;det++)
      {
	Int_t nRings = (det==1 ? 1 : 2);
	for(Int_t ring = 0;ring<nRings;ring++)
	  {
	    Char_t ringChar = (ring == 0 ? 'I' : 'O');
	    for(Int_t v=0; v<fNvtxBins;v++)
	      {
		
		TH2F* hHits          = (TH2F*)listOfHits->FindObject(Form("hHits_FMD%d%c_vtx%d", det,ringChar,v));
		fListOfHits.Add(hHits);
	      }
	  }
      }
  for(Int_t iring = 0; iring<2;iring++) {
    Char_t ringChar = (iring == 0 ? 'I' : 'O');
    for(Int_t v=0; v<fNvtxBins;v++) {
      
      TH2F* hPrimary       = (TH2F*)listOfPrim->FindObject( Form("hPrimary_FMD_%c_vtx%d",ringChar,v));
      fListOfPrimaries.Add(hPrimary);
      
    }
  }
  GenerateCorrection();
  
  TFile fout("backgroundFromFile.root","recreate");
  fListOfHits.Write();
  fListOfPrimaries.Write();
  fListOfCorrection.Write();
  fVertexBins.Write();
  fout.Close();
  
  if(storeInOCDB) {
    TFile fcalib("$ALICE_ROOT/FMD/Correction/Background/background.root","RECREATE");
    fBackground->Write(AliFMDAnaParameters::fkBackgroundID);
    fcalib.Close();
    /*  AliCDBManager* cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    AliCDBId      id(AliFMDAnaParameters::GetBackgroundPath(),runNo,999999999);
    
    AliCDBMetaData* meta = new AliCDBMetaData;				
    meta->SetResponsible(gSystem->GetUserInfo()->fRealName.Data());	
    meta->SetAliRootVersion(gROOT->GetVersion());			
    meta->SetBeamPeriod(1);						
    meta->SetComment("Background Correction for FMD");
    meta->SetProperty("key1", fBackground );
    cdb->Put(fBackground, id, meta);
    */
  }
  
}
//_____________________________________________________________________
//
// EOF
//
