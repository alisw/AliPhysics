//SSD dEdx gain calibration task
//to run in QA train
//autor Marek Chojnacki
//Marek.Chojnacki@cern.ch
// last line
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliITSCalibrationSSD.h"

#include "AliAnalysisTaskdEdxSSDQA.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"


#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "Riostream.h"

#include "AliESDfriend.h"
#include "AliESDfriendTrack.h"
#include "AliTrackPointArray.h"
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliITSgeomTGeo.h"
#include "TMath.h"
#include "TGeoMatrix.h"

using namespace std;
ClassImp(AliAnalysisTaskdEdxSSDQA)

//________________________________________________________________________
AliAnalysisTaskdEdxSSDQA::AliAnalysisTaskdEdxSSDQA(const char *name) 
: AliAnalysisTaskSE(name), 
  fHistBadMods(0),
  fHistBadPstr(0),
  fHistBadNstr(0),
  fHist1(0),
  fHist2(0),
  fHist3(0),
  fHist4(0),
  fHist5(0),
  fHist6(0),
  fHist7(0),
  fHist8(0),
  fHist9(0),
  fHist10(0),
  fHist1sa(0),
  fHist2sa(0),
  fHist3sa(0),
  fHist4sa(0),
  fHist5sa(0),
  fHist6sa(0),
  fHist7sa(0),
  fHist8sa(0),
  fHist9sa(0),
  fHist10sa(0),
  fListOfHistos(0),
  fPcut(0.0),
  fdothecorrection(0),
  fInitCalib(kFALSE),
  fSSDCalibration(0x0)
{
  // Constructor
  // Define input and output slots here
  // Input slot #1 works with a TChain
  // fcorrections=new Float_t[1698*12];
  for(int i=0;i<1698*12;i++)
    fcorrections[i]=1.0;
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskdEdxSSDQA::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  	
  fListOfHistos = new TList();
  fListOfHistos->SetOwner();

  fHistBadMods=new TH1F("BadModules","",1698,-0.5,1697.5);
  fListOfHistos->Add(fHistBadMods);

  fHistBadPstr=new TH1F("Bad-p-strips","",1698,-0.5,1697.5);
  fListOfHistos->Add(fHistBadPstr);

  fHistBadNstr=new TH1F("Bad-n-strips","",1698,-0.5,1697.5);
  fListOfHistos->Add(fHistBadNstr);

  fHist1 =new TH2F("QAChargeRatio","QAChargeRatio;Module;CR",1698,-0.5,1697.5,80,-1.0,1.0);
  fListOfHistos->Add(fHist1);
  fHist2=new TH2F("QACharge","QACharge;Module;Q",1698,-0.5,1697.5,150,0,300);
  fListOfHistos->Add(fHist2);
  fHist3=new TH3F("ChargeRatiovCharge","CRvQ;module;Q;CR",1698,-0.5,1697.5,80,0,1600,50,-1.0,1.0);
  fListOfHistos->Add(fHist3);
  fHist4=new TH2F("QAChargedinchips","chargedinchips;chip;Q",1698*12,-0.5,1698*12-0.5,150,0,300);
  fListOfHistos->Add(fHist4);

  fHist1sa =new TH2F("QAChargeRatioSA","QAChargeRatio;Module;CR",1698,-0.5,1697.5,80,-1.0,1.0);
  fListOfHistos->Add(fHist1sa);
  fHist2sa=new TH2F("QAChargeSA","QACharge;Module;Q",1698,-0.5,1697.5,150,0,300);
  fListOfHistos->Add(fHist2sa);
  fHist3sa=new TH3F("ChargeRatiovChargeSA","CRvQ;module;Q;CR",1698,-0.5,1697.5,80,0,1600,50,-1.0,1.0);
  fListOfHistos->Add(fHist3sa);
  fHist4sa=new TH2F("QAChargedinchipsSA","chargedinchips;chip;Q",1698*12,-0.5,1698*12-0.5,150,0,300);
  fListOfHistos->Add(fHist4sa);
  if(fdothecorrection)
    {
      fHist5=new TH2F("QAChargeCorrected","chargedCorrected;Module;QCorrected",1698,-0.5,1697.5,150,0,300);
      fListOfHistos->Add(fHist5);
      fHist5sa=new TH2F("QAChargeCorrectedSA","chargedCorrected;Module;QCorrected",1698,-0.5,1697.5,150,0,300);
      fListOfHistos->Add(fHist5sa);
    }
  fHist6=new TH2F("QNvQP",";QN;QP",160,0,1600,1600,0,1600);
  fListOfHistos->Add(fHist6);
  fHist6sa=new TH2F("QNvQPSA",";QN;QP",160,0,1600,1600,0,1600);
  fListOfHistos->Add(fHist6sa);
	
  fHist7=new TH2F("ZetaPhiLay5","Layer 5 - global tracks ; z (cm) ; #varphi (rad)",100,-50,50,100,0.,2.*TMath::Pi());
  fListOfHistos->Add(fHist7);
  fHist7sa=new TH2F("ZetaPhiLay5sa","Layer 5 - ITSsa tracks ; z (cm) ; #varphi (rad)",100,-50,50,100,0.,2.*TMath::Pi());
  fListOfHistos->Add(fHist7sa);

  fHist8=new TH2F("ZetaPhiLay6","Layer 6 - global tracks ; z (cm) ; #varphi (rad)",100,-50,50,100,0.,2.*TMath::Pi());
  fListOfHistos->Add(fHist8);
  fHist8sa=new TH2F("ZetaPhiLay6sa","Layer 6 - ITSsa tracks ; z (cm) ; #varphi (rad)",100,-50,50,100,0.,2.*TMath::Pi());
  fListOfHistos->Add(fHist8sa);

  fHist9=new TH2F("LadModLay5","Layer 5 - global tracks ; Detector ; Ladder",22,0.5,22.5,34,0.5,34.5);
  fListOfHistos->Add(fHist9);
  fHist9sa=new TH2F("LadModLay5sa","Layer 5 - ITSsa tracks ; Detector ; Ladder",22,0.5,22.5,34,0.5,34.5);
  fListOfHistos->Add(fHist9sa);

  fHist10=new TH2F("LadModLay6","Layer 6 - global tracks ; Detector ; Ladder",25,0.5,25.5,38,0.5,38.5);
  fListOfHistos->Add(fHist10);
  fHist10sa=new TH2F("LadModLay6sa","Layer 6 - ITSsa tracks ;  Detector ; Ladder",25,0.5,25.5,38,0.5,38.5);
  fListOfHistos->Add(fHist10sa);


  PostData(1,  fListOfHistos);
}
//________________________________________________________________________
void AliAnalysisTaskdEdxSSDQA::LocalInit() 
{
	
  Printf("end of LocalInit");
}
//______________________________________________________________________
AliAnalysisTaskdEdxSSDQA::~AliAnalysisTaskdEdxSSDQA()
{
  // desctructor: no need to delete histos because the TList is owner
  delete fListOfHistos;
}

//________________________________________________________________________

void AliAnalysisTaskdEdxSSDQA::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  AliESDEvent* esd = dynamic_cast<AliESDEvent*> (InputEvent());
  if(!esd) 
    {
      Printf("ERROR: Input ESD Event not available");
      PostData(1,  fListOfHistos);
      return;
    }
    
  Int_t nTracks=esd->GetNumberOfTracks();

  if(!ESDfriend())
    {
      Printf("problem with friend");
      PostData(1,  fListOfHistos);
      return;
    }

  if (!fInitCalib) {
    AliCDBManager* man = AliCDBManager::Instance();
    if (!man) {
       AliFatal("CDB not set but needed by AliAnalysisTaskSDDRP");
       return;
    }   
    fSSDCalibration=new AliITSCalibrationSSD();
    Bool_t okCalib=kTRUE;
    AliCDBEntry* e1=(AliCDBEntry*)man->Get("ITS/Calib/BadChannelsSSD");
    if(e1){
      e1->PrintId();
      e1->PrintMetaData();
      AliITSBadChannelsSSDv2* badChannelsSSD=(AliITSBadChannelsSSDv2*)e1->GetObject();
      fSSDCalibration->SetBadChannels(badChannelsSSD);
    }else{
      printf("BadChannelsSSD object not found\n");
      okCalib=kFALSE;
    }
    AliCDBEntry* e2=(AliCDBEntry*)man->Get("ITS/Calib/NoiseSSD");
    if(e2){
      e2->PrintId();
      e2->PrintMetaData();
      AliITSNoiseSSDv2* noiseSSD=(AliITSNoiseSSDv2*)e2->GetObject();
      fSSDCalibration->SetNoise(noiseSSD);
    }else{
      printf("NoiseSSD object not found\n");
      okCalib=kFALSE;
    }
    AliCDBEntry* e3=(AliCDBEntry*)man->Get("ITS/Calib/GainSSD");
    if(e3){
      e3->PrintId();
      e3->PrintMetaData();
      AliITSGainSSDv2* gainSSD=(AliITSGainSSDv2*)e3->GetObject();
      fSSDCalibration->SetGain(gainSSD);
    }else{
      printf("GainSSD object not found\n");
      okCalib=kFALSE;
    }
    if(okCalib){
      for(Int_t imod=0; imod<1698; imod++){
	fSSDCalibration->SetModule(imod);
	if(fSSDCalibration->IsBad()) fHistBadMods->SetBinContent(imod+1,0.);
	else fHistBadMods->SetBinContent(imod+1,1.);
	Int_t countPbad=0, countNbad=0;
	for(Int_t ib=0; ib<768; ib++) {
	  if(fSSDCalibration->IsPChannelBad(ib)) countPbad++;
	  if(fSSDCalibration->IsNChannelBad(ib)) countNbad++;
	}
	fHistBadPstr->SetBinContent(imod+1,countPbad);
	fHistBadNstr->SetBinContent(imod+1,countNbad);
      }
    }
    fInitCalib = kTRUE;
  }  
    

  //Printf("Event nTracks %d", nTracks);	
  AliTrackPointArray*  trackar=0x0;
  Bool_t l5;
  Bool_t l6;
  Int_t npoints;
  AliTrackPoint point;
  Int_t nlayer=0;
  Int_t id=0;
  Int_t iLayer=0;
  Int_t iLadder=0;
  Int_t iDetec=0;
  Float_t chargeratio=0.0;
  Float_t charge=0.0;
  for(int itr=0;itr<nTracks;itr++)
    {
      AliESDtrack* track= esd->GetTrack(itr);

      if(TMath::Abs(track->Eta())>0.9)
	continue;
      if(track->GetP()>10.0)
	continue;	
      if (track->IsOn(AliESDtrack::kITSrefit)){
	Bool_t goodGlobal=kFALSE;
	Bool_t goodSA=kFALSE;
	if(track->IsOn(AliESDtrack::kTPCrefit)) goodGlobal=kTRUE;
	if(!(track->IsOn(AliESDtrack::kTPCin))) goodSA=kTRUE;

	l5=track->HasPointOnITSLayer(4);
	l6=track->HasPointOnITSLayer(5);
	if (!(l5||l6))//only tracks with SSD point
	  continue;
	Double_t tmpQESD[4]={-1.0,-1.0,-1.0,-1.0};
	track->GetITSdEdxSamples(tmpQESD);
	trackar=(AliTrackPointArray*)track->GetTrackPointArray(); 
	if(!trackar)
	  continue;
	npoints=trackar->GetNPoints();			
	for(int itnp=0;itnp<npoints;itnp++)		
	  {									
	    if(trackar->GetPoint(point,itnp))
	      {
		nlayer=AliGeomManager::VolUIDToLayer(point.GetVolumeID(), id);//layer number
	      }
	    else
	      continue;		
	    if(nlayer==5||nlayer==6)
	      {
		TGeoHMatrix* geomatrix=AliGeomManager::GetMatrix(point.GetVolumeID());
		//geomatrix->Print();
		if(point.GetCharge()>0.0&&point.IsExtra()==kFALSE)//do not use additional clusters
		  {		
		    Double_t local[3]={0.0,0.0,0.0};
		    Double_t global[3]={0.0,0.0,0.0};
		    global[0]=point.GetX();
		    global[1]=point.GetY();
		    global[2]=point.GetZ();
		    Double_t phi=TMath::ATan2(-global[1],-global[0])+TMath::Pi();
		    Int_t modId=id+AliITSgeomTGeo::GetModuleIndex(nlayer,1,1);
		    AliITSgeomTGeo::GetModuleId(modId,iLayer,iLadder,iDetec);
		    geomatrix->MasterToLocal(global,local);
		    chargeratio=point.GetChargeRatio();	
		    charge=point.GetCharge();
		    Float_t fQNnotcorr=charge*(1.0+chargeratio);
		    Float_t fQPnotcorr=charge*(1.0-chargeratio);						 
		    if(goodGlobal) fHist6->Fill(fQNnotcorr,fQPnotcorr);
		    if(goodSA) fHist6sa->Fill(fQNnotcorr,fQPnotcorr);
		    // cout<<point.GetCharge()<<" "<<point.GetChargeRatio()<<" "<<nlayer<<" "<<id<<endl;
		    if(nlayer==5&&tmpQESD[2]>0.0)
		      {
			if(goodGlobal){ 
			  fHist1->Fill(id,chargeratio);
			  fHist3->Fill(id,charge,chargeratio);
			  fHist7->Fill(global[2],phi);
			  fHist9->Fill(iDetec,iLadder);
			}
			if(goodSA){ 
			  fHist1sa->Fill(id,chargeratio);
			  fHist3sa->Fill(id,charge,chargeratio);
			  fHist7sa->Fill(global[2],phi);
			  fHist9sa->Fill(iDetec,iLadder);
			}
							
			if(track->GetP()>fPcut)
			  {
			    Float_t fQP=tmpQESD[2]*(1.0-chargeratio);
			    Float_t fQN=tmpQESD[2]*(1.0+chargeratio);
			    Int_t fPchip=Pstrip5(10000.0*local[0],10000.0*local[2])/128;
			    Int_t fNchip=Nstrip5(10000.0*local[0],10000.0*local[2])/128;
								
			    if(goodGlobal){
			      fHist2->Fill(id,tmpQESD[2]);
			      fHist4->Fill(id*12+fPchip,fQP);	
			      fHist4->Fill(id*12+6+fNchip,fQN);
			    }
			    if(goodSA){
			      fHist2sa->Fill(id,tmpQESD[2]);
			      fHist4sa->Fill(id*12+fPchip,fQP);	
			      fHist4sa->Fill(id*12+6+fNchip,fQN);
			    }
			    if(fdothecorrection){
			      if(goodGlobal) fHist5->Fill(id,0.5*(fQP*fcorrections[12*id+fPchip]+fQN*fcorrections[12*id+6+fNchip]));	
			      if(goodSA) fHist5sa->Fill(id,0.5*(fQP*fcorrections[12*id+fPchip]+fQN*fcorrections[12*id+6+fNchip]));	
			    }
			  }	
		      }	
		    if(nlayer==6&&tmpQESD[3]>0.0)
		      {
			if(goodGlobal){ 
			  fHist1->Fill(id+748,chargeratio);
			  fHist3->Fill(id+748,charge,chargeratio);	
			  fHist8->Fill(global[2],phi);
			  fHist10->Fill(iDetec,iLadder);
			}
			if(goodSA){ 
			  fHist1sa->Fill(id+748,chargeratio);
			  fHist3sa->Fill(id+748,charge,chargeratio);	
			  fHist8sa->Fill(global[2],phi);
			  fHist10sa->Fill(iDetec,iLadder);
			}
			if(track->GetP()>fPcut)
			  {
							
			    Float_t fQP=tmpQESD[3]*(1.0-chargeratio);
			    Float_t fQN=tmpQESD[3]*(1.0+chargeratio);
			    Int_t fPchip=Pstrip6(10000.0*local[0],10000.0*local[2])/128;
			    Int_t fNchip=Nstrip6(10000.0*local[0],10000.0*local[2])/128;

			    if(goodGlobal){							
			      fHist2->Fill(id+748,tmpQESD[3]);
			      fHist4->Fill((id+748)*12+fPchip,fQP);	
			      fHist4->Fill((id+748)*12+6+fNchip,fQN);
			    }
			    if(goodSA){
			      fHist2sa->Fill(id+748,tmpQESD[3]);
			      fHist4sa->Fill((id+748)*12+fPchip,fQP);	
			      fHist4sa->Fill((id+748)*12+6+fNchip,fQN);
			    }
			    if(fdothecorrection){
			      if(goodGlobal) fHist5->Fill((id+748),0.5*(fQP*fcorrections[12*(id+748)+fPchip]+fQN*fcorrections[12*(id+748)+6+fNchip]));			
			      if(goodSA) fHist5sa->Fill((id+748),0.5*(fQP*fcorrections[12*(id+748)+fPchip]+fQN*fcorrections[12*(id+748)+6+fNchip]));			
			    }
			  }	
		      }
							
		  }	
	      }
	  }	
      }
    }
  //Printf("Event nTracks end %d", nTracks);
  // Post output data.
  PostData(1,  fListOfHistos);
}      

//________________________________________________________________________
void AliAnalysisTaskdEdxSSDQA::Terminate(Option_t *) 
{
  //terminate
  Printf("end of Terminate");
}
//_____________________________________________________________________________
Int_t  AliAnalysisTaskdEdxSSDQA::Pstrip5(Float_t x,Float_t z) const
{
  // P strip for layer 5 from local cooridnates 
  Float_t value;
  x=x+36432.5-223.9;
  z=z-(76000.0/7.0);
  value=-(3.0/38000.0)*z+x/95.0;
  return TMath::Nint(value);
}
//___________________________________________________________________________
Int_t  AliAnalysisTaskdEdxSSDQA::Nstrip5(Float_t x,Float_t z)  const
{
  // N strip for layer 5 from local cooridnates 
  Float_t value;
  x=x+36432.5-223.9;
  z=z-(76000.0/7.0);
  value=x/95.0+(11.0/38000.0)*z;
  return TMath::Nint(value);
}
//___________________________________________________________________________
Int_t  AliAnalysisTaskdEdxSSDQA::Pstrip6(Float_t x,Float_t z)  const
{
  // P strip for layer 6 from local cooridnates 
  Float_t value;
  x=x-36432.5-223.9;
  z=z+(76000.0/7.0);
  value=-(3.0/38000.0)*z-x/95.0;
  return TMath::Nint(value);
}
//__________________________________________________________________________________
Int_t  AliAnalysisTaskdEdxSSDQA::Nstrip6(Float_t x,Float_t z)  const
{
  // N strip for layer 6 from local cooridnates 
  Float_t value;
  x=x-36432.5-223.9;
  z=z+(76000.0/7.0);
  value=x*(-1.0/95.0)+(11.0/38000.0)*z;
  return TMath::Nint(value);
}
//___________________________________________________________________________________
void  AliAnalysisTaskdEdxSSDQA::SetDoChipCorretions(const char* filename)
{
  //Upload corrections for each chip only for test 
  cout<<filename<<endl;
  ifstream infile(filename);
  if(!infile.good())
    return;
  Float_t value=0.0;
  fdothecorrection=1;
  for (int i=0;i<1698;i++)
    {
      for (int j=0;j<6;j++)
	{
	  infile>>value;
	  cout<<value<<" ";
	  if(value>0.0)
	    fcorrections[i*12+j]=84.0/value;
	  else
	    fcorrections[i*12+j]=1.0;	
	  infile>>value;
	  infile>>value;
	  cout<<value<<" ";
	  if(value>0.0)
	    fcorrections[i*12+j+6]=84.0/value;
	  else
	    fcorrections[i*12+j+6]=1.0;
	  infile>>value;
	}
      cout<<endl;		
    }	
}

