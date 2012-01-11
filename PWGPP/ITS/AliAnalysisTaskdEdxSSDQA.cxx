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
#include "TMath.h"
#include "TGeoMatrix.h"

using namespace std;
ClassImp(AliAnalysisTaskdEdxSSDQA)

//________________________________________________________________________
AliAnalysisTaskdEdxSSDQA::AliAnalysisTaskdEdxSSDQA(const char *name) 
: AliAnalysisTaskSE(name), fHist1(0),fHist2(0),fHist3(0),fHist4(0),fHist5(0),fHist6(0),fListOfHistos(0),fPcut(0.0),fdothecorrection(0)
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
	fHist1 =new TH2F("QAChargeRatio","QAChargeRatio;Module;CR",1698,-0.5,1697.5,80,-1.0,1.0);
	fListOfHistos->Add(fHist1);
	fHist2=new TH2F("QACharge","QACharge;Module;Q",1698,-0.5,1697.5,150,0,300);
	fListOfHistos->Add(fHist2);
	fHist3=new TH3F("ChargeRatiovCharge","CRvQ;module;Q;CR",1698,-0.5,1697.5,80,0,1600,50,-1.0,1.0);
	fListOfHistos->Add(fHist3);
	fHist4=new TH2F("QAChargedinchips","chargedinchips;chip;Q",1698*12,-0.5,1698*12-0.5,150,0,300);
	fListOfHistos->Add(fHist4);
	if(fdothecorrection)
	{
		fHist5=new TH2F("QAChargeCorrected","chargedCorrected;Module;QCorrected",1698,-0.5,1697.5,150,0,300);
		fListOfHistos->Add(fHist5);
	}
	fHist6=new TH2F("QNvQP",";QN;QP",160,0,1600,1600,0,1600);
	fListOfHistos->Add(fHist6);
	
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
	//Printf("Event nTracks %d", nTracks);	
	AliTrackPointArray*  trackar=0x0;
	Bool_t l5;
	Bool_t l6;
	Int_t npoints;
	AliTrackPoint point;
	Int_t nlayer=0;
	Int_t id=0;
	Float_t chargeratio=0.0;
	Float_t charge=0.0;
	for(int itr=0;itr<nTracks;itr++)
	{
		AliESDtrack* track= esd->GetTrack(itr);

		if(TMath::Abs(track->Eta())>0.9)
			continue;
		if(track->GetP()>10.0)
			continue;	
		if (track->IsOn(AliESDtrack::kITSrefit)&&track->IsOn(AliESDtrack::kTPCrefit))
		{
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
						geomatrix->MasterToLocal(global,local);
						 chargeratio=point.GetChargeRatio();	
						 charge=point.GetCharge();
						 Float_t fQNnotcorr=charge*(1.0+chargeratio);
						 Float_t fQPnotcorr=charge*(1.0-chargeratio);						 
						 fHist6->Fill(fQNnotcorr,fQPnotcorr);
						// cout<<point.GetCharge()<<" "<<point.GetChargeRatio()<<" "<<nlayer<<" "<<id<<endl;
						if(nlayer==5&&tmpQESD[2]>0.0)
						{
							fHist1->Fill(id,chargeratio);
							if(track->GetP()>fPcut)
							{
								Float_t fQP=tmpQESD[2]*(1.0-chargeratio);
								Float_t fQN=tmpQESD[2]*(1.0+chargeratio);
								Int_t fPchip=Pstrip5(10000.0*local[0],10000.0*local[2])/128;
								Int_t fNchip=Nstrip5(10000.0*local[0],10000.0*local[2])/128;
								
								fHist2->Fill(id,tmpQESD[2]);
								fHist4->Fill(id*12+fPchip,fQP);	
								fHist4->Fill(id*12+6+fNchip,fQN);
								if(fdothecorrection)
									fHist5->Fill(id,0.5*(fQP*fcorrections[12*id+fPchip]+fQN*fcorrections[12*id+6+fNchip]));	
							}	
							fHist3->Fill(id,charge,chargeratio);
						}	
						if(nlayer==6&&tmpQESD[3]>0.0)
						{
							fHist1->Fill(id+748,chargeratio);
							if(track->GetP()>fPcut)
							{
							
								Float_t fQP=tmpQESD[3]*(1.0-chargeratio);
								Float_t fQN=tmpQESD[3]*(1.0+chargeratio);
								Int_t fPchip=Pstrip6(10000.0*local[0],10000.0*local[2])/128;
								Int_t fNchip=Nstrip6(10000.0*local[0],10000.0*local[2])/128;
							
								fHist2->Fill(id+748,tmpQESD[3]);
								fHist4->Fill((id+748)*12+fPchip,fQP);	
								fHist4->Fill((id+748)*12+6+fNchip,fQN);
								if(fdothecorrection)
									fHist5->Fill((id+748),0.5*(fQP*fcorrections[12*(id+748)+fPchip]+fQN*fcorrections[12*(id+748)+6+fNchip]));			
								
							}	
							fHist3->Fill(id+748,charge,chargeratio);	
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

