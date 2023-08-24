///
/// \file AliFemtoCutMonitorDphiDeta1.cxx
///

#include "AliFemtoCutMonitorDphiDeta1.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoEvent.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>

//fPtmax=0.0;
//fPhimax=0.0;
//fEtamax=0.0;


AliFemtoCutMonitorDphiDeta1::AliFemtoCutMonitorDphiDeta1():
  fDPhiDEta(NULL),
  fPtmax(0),
  fPhimax(0),
  fEtamax(0)
{
  // Default constructor
  fDPhiDEta = new TH2D("fDPhiDEta","fDPhiDEta",
			100, -3.0, 3.0,
                        100, -6.0, 6.0);
}

AliFemtoCutMonitorDphiDeta1::AliFemtoCutMonitorDphiDeta1(const char *aName, int nBins, double multMax):
  AliFemtoCutMonitor()
{
  TString name(aName);

  // Normal constructor

  fDPhiDEta = new TH2D("fDPhiDEta"+name, "fDPhiDEta",
			100, -3.0, 3.0,
                        100, -6.0, 6.0);
 

}

AliFemtoCutMonitorDphiDeta1::AliFemtoCutMonitorDphiDeta1(const AliFemtoCutMonitorDphiDeta1 &aCut):
  AliFemtoCutMonitor(aCut),
  fDPhiDEta(NULL)
{
  // copy constructor
  fDPhiDEta = new TH2D(*aCut.fDPhiDEta);
  

}

AliFemtoCutMonitorDphiDeta1::~AliFemtoCutMonitorDphiDeta1()
{
  // Destructor
  delete fDPhiDEta ;
}

AliFemtoCutMonitorDphiDeta1& AliFemtoCutMonitorDphiDeta1::operator=(const AliFemtoCutMonitorDphiDeta1& aCut)
{
  // assignment operator
  if (this == &aCut)
    return *this;


  if (fDPhiDEta) delete fDPhiDEta ;
  fDPhiDEta  = new TH2D(*aCut.fDPhiDEta);
 

  return *this;
}

AliFemtoString AliFemtoCutMonitorDphiDeta1::Report()
{
  // Prepare report from the execution
  TString report = "*** AliFemtoCutMonitorDphiDeta1 report";
  return AliFemtoString((const char *)report);
}

void AliFemtoCutMonitorDphiDeta1::Fill(const AliFemtoEvent* aEvent)
{
  // Fill in the monitor histograms with the values from the current track
  

 int mult = (int) aEvent->UncorrectedNumberOfPrimaries();


  Int_t ParticleNumber = 0;

   AliFemtoTrackCollection * tracks = aEvent->TrackCollection(); 
   
   Double_t Ptmax=0; 
   Double_t Phimax=0; 
   Double_t Etamax=0; 
   Double_t DEta=0;
   Double_t DPhi=0;
    
  for (AliFemtoTrackIterator iter=tracks->begin();iter!=tracks->end();iter++){

    Double_t NewPhi = (*iter)->P().Phi();
    Double_t NewPt =  (*iter)->Pt();
    Double_t NewEta = (*iter)->P().PseudoRapidity();

    if(NewPt > Ptmax){
      Ptmax = NewPt;
      Phimax = NewPhi;
      Etamax = NewEta;
    }
    ParticleNumber++;

  }  // end of track loop


  fPtmax = Ptmax;
  fPhimax = Phimax;
  fEtamax = Etamax;

 if(Ptmax>0.5){
   for (AliFemtoTrackIterator iter=tracks->begin();iter!=tracks->end();iter++){

    Double_t NewPhi = (*iter)->P().Phi();
    Double_t NewPt =  (*iter)->Pt();
    Double_t NewEta = (*iter)->P().PseudoRapidity();
    
    if(TMath::Abs(Ptmax-NewPt)>0.001){

    DEta = NewEta - Etamax;
    DPhi = NewPhi - Phimax;

    ParticleNumber++;

    fDPhiDEta->Fill(DEta,DPhi);

     }

    }

  }  // end of track loop

}

void AliFemtoCutMonitorDphiDeta1::Write()
{
  // Write out the relevant histograms
  fDPhiDEta->Write();
}

TList *AliFemtoCutMonitorDphiDeta1::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fDPhiDEta);

  return tOutputList;
}

