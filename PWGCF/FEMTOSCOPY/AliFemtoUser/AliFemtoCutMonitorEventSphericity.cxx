///
/// \file AliFemtoCutMonitorEventSphericity.cxx
///

#include "AliFemtoCutMonitorEventSphericity.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoEvent.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>

AliFemtoCutMonitorEventSphericity::AliFemtoCutMonitorEventSphericity():
  fEvSpher(NULL),
  fEvSpherMult(NULL),
  fMultSumPt(NULL)
{
  // Default constructor
  fEvSpher = new TH1D("EvSpher", "EvSpher", 100, 0., 1.0);
  fEvSpherMult = new TH2D("EvSpherMult","EvSpher vs Mult",
			100, 0., 1.0,
                        5001, -0.5, 5000.5);
  fMultSumPt = new TH2D("EvMultSumPt",
                        "Event Multiplicity vs Total pT",
                        5001, -0.5, 5000.5,
                        1000, 0.0, 100.0);
}

AliFemtoCutMonitorEventSphericity::AliFemtoCutMonitorEventSphericity(const char *aName, int nBins, double multMax):
  AliFemtoCutMonitor(),
  fEvSpher(NULL),
  fEvSpherMult(NULL),
  fMultSumPt(NULL)
{
  TString name(aName);

  // Normal constructor
 
  
  fEvSpher = new TH1D("EvSpher" + name,
                     "Event Spher",
                      100, 0., 1.0);

  fEvSpherMult = new TH2D("EvSpherMult"+ name, "EvSpher vs Mult",
			100, 0., 1.0,
                        5001, -0.5, 5000.5);


  fMultSumPt = new TH2D("EvMultTotPt" + name,
                        "Event Multiplicity vs Total pT",
                        501, -0.5, 500.5,
                        1000, 0.0, 100.0);


}

AliFemtoCutMonitorEventSphericity::AliFemtoCutMonitorEventSphericity(const AliFemtoCutMonitorEventSphericity &aCut):
  AliFemtoCutMonitor(aCut),
  fEvSpher(NULL),
  fEvSpherMult(NULL),
  fMultSumPt(NULL)
{
  // copy constructor
  fEvSpher = new TH1D(*aCut.fEvSpher);
  fEvSpherMult = new TH2D(*aCut.fEvSpherMult);
  fMultSumPt = new TH2D(*aCut.fMultSumPt);

}

AliFemtoCutMonitorEventSphericity::~AliFemtoCutMonitorEventSphericity()
{
  // Destructor
  delete fEvSpher;
  delete fEvSpherMult;
  delete fMultSumPt;
}

AliFemtoCutMonitorEventSphericity& AliFemtoCutMonitorEventSphericity::operator=(const AliFemtoCutMonitorEventSphericity& aCut)
{
  // assignment operator
  if (this == &aCut)
    return *this;

  if (fEvSpher) delete fEvSpher;
  fEvSpher = new TH1D(*aCut.fEvSpher);

  if (fEvSpherMult) delete fEvSpherMult;
  fEvSpherMult = new TH2D(*aCut.fEvSpherMult);
 
  if (fMultSumPt) delete fMultSumPt;
  fMultSumPt = new TH2D(*aCut.fMultSumPt);


  return *this;
}

AliFemtoString AliFemtoCutMonitorEventSphericity::Report()
{
  // Prepare report from the execution
  TString report = "*** AliFemtoCutMonitorEventSphericity report";
  return AliFemtoString((const char *)report);
}

void AliFemtoCutMonitorEventSphericity::Fill(const AliFemtoEvent* aEvent)
{
  // Fill in the monitor histograms with the values from the current track
  

 int mult = (int) aEvent->UncorrectedNumberOfPrimaries();


  Int_t ParticleNumber = 0;
  Double_t SumPt = 0;
  Double_t S00=0; 
  Double_t S11=0;
  Double_t S10=0;
  Double_t Lambda1 = 0;
  Double_t Lambda2 = 0;
  Double_t St = 0;

   AliFemtoTrackCollection * tracks = aEvent->TrackCollection(); 
   
   
  for (AliFemtoTrackIterator iter=tracks->begin();iter!=tracks->end();iter++){
  
  
    Double_t NewPhi = (*iter)->P().Phi();
    Double_t NewPt =  (*iter)->Pt();
    Double_t NewEta = (*iter)->P().PseudoRapidity();
   
    
    if(TMath::Abs(NewEta)>0.8 || NewPt<0.5){continue;}
    
    Double_t Px;
    Double_t Py;
    
    Px= NewPt * TMath::Cos(NewPhi);
    Py= NewPt * TMath::Sin(NewPhi);
    
    S00 = S00 + Px*Px/(NewPt);  // matrix elements of the transverse shpericity matrix S(i,j)
    S11 = S11 + Py*Py/(NewPt);  // i,j /in [0,1]
    S10 = S10 + Px*Py/(NewPt);
    SumPt = SumPt + NewPt;
    ParticleNumber++;
    
  }  	// end of track loop

    if(SumPt==0){St=0; goto M1;}
      
  S00 = S00/SumPt; // normalize
  S11 = S11/SumPt;
  S10 = S10/SumPt;
  
  Lambda1 = (S00 + S11 + TMath::Sqrt((S00+S11)*(S00+S11)-4.0*(S00*S11-S10*S10)))/2.0;
  Lambda2 = (S00 + S11 - TMath::Sqrt((S00+S11)*(S00+S11)-4.0*(S00*S11-S10*S10)))/2.0;
  
     if(Lambda1+Lambda2!=0 && ParticleNumber>2)
	{
		St = 2*Lambda2/(Lambda1+Lambda2);
	}
     else{St=0; goto M1;};
  

 M1: // cout<<"St  = "<< St << " " <<mult<<" "<<SumPt<<endl;

  fEvSpher->Fill(St);
  fEvSpherMult->Fill(St,mult);
  fMultSumPt->Fill(mult,SumPt);



}

void AliFemtoCutMonitorEventSphericity::Write()
{
  // Write out the relevant histograms
  fEvSpher->Write();
  fEvSpherMult->Write();
  fMultSumPt->Write();

}

TList *AliFemtoCutMonitorEventSphericity::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fEvSpher);
  tOutputList->Add(fEvSpherMult);
  tOutputList->Add(fMultSumPt);

  return tOutputList;
}

//void AliFemtoCutMonitorEventSphericity::SetReadMC(Bool_t mc)
//{
//  freadMC = mc;
//}

//void AliFemtoCutMonitorEventSphericity::AdditionalMultHistsOn(Bool_t addhists)
//{
//  faddhists = addhists;
//}
