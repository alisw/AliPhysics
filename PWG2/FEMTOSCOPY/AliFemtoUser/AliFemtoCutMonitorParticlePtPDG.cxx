////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCutMonitorParticlePtPDG - the cut monitor for particles to study    //
// the difference between reconstructed and true momentum                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitorParticlePtPDG.h"
#include "AliFemtoModelHiddenInfo.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TMath.h>

AliFemtoCutMonitorParticlePtPDG::AliFemtoCutMonitorParticlePtPDG():
  fPtPDG(0),ftpcHist(0),fPtGoodPi(0),fPtFakePi(0),fPtGoodK(0),fPtFakeK(0),
  fPtGoodP(0),fPtFakeP(0),fPtRPi(0),fPtRK(0),fPtRP(0),
  fPtContP(0),
  fPtContPi(0),
  fPtContMup(0),
  fPtContElp(0),
  fMass(0.13957)
{
  // Default constructor
}

AliFemtoCutMonitorParticlePtPDG::AliFemtoCutMonitorParticlePtPDG(const char *aName, float aMass):
  AliFemtoCutMonitor(),
  fPtPDG(0),ftpcHist(0),fPtGoodPi(0),fPtFakePi(0),fPtGoodK(0),fPtFakeK(0),
  fPtGoodP(0),fPtFakeP(0),fPtRPi(0),fPtRK(0),fPtRP(0),
  fPtContP(0),
  fPtContPi(0),
  fPtContMup(0),
  fPtContElp(0),
  fMass(aMass)
{
  // Normal constructor
  char name[200];
  snprintf(name, 200, "PtPDG%s", aName);
  fPtPDG = new TH2D(name, "PDG vs Pt", 10, 0.0, 5.0, 100, 0.1, 2.0);
  snprintf(name, 200, "tpcHist%s", aName);
  ftpcHist=new TH2D(name,"TPC dE/dX vs momentum",100,0.1,2.7,100,0.,6.);
  snprintf(name, 200, "PtGoodPi%s", aName);
  fPtGoodPi = new TH1D(name, "good pions Pt",                    100, 0.1, 2.0);
  snprintf(name, 200, "PtFakePi%s", aName);
  fPtFakePi = new TH1D(name, "fake pions Pt",              100, 0.1, 2.0);
  snprintf(name, 200, "PtRPi%s", aName);
  fPtRPi = new TH1D(name, "right pdg pions Pt",               100, 0.1, 2.0);
  snprintf(name, 200, "PtGoodK%s", aName);
  fPtGoodK = new TH1D(name, "good kaons Pt",                     100, 0.1, 2.0);
  snprintf(name, 200, "PtFakeK%s", aName);
  fPtFakeK = new TH1D(name, "fake kaons Pt",                100, 0.1, 2.0);
  snprintf(name, 200, "PtRK%s", aName);
  fPtRK = new TH1D(name, "right pdg kaons Pt",                 100, 0.1, 2.0);  
   snprintf(name, 200, "PtGoodP%s", aName);
  fPtGoodP = new TH1D(name, "good protons Pt",                     100, 0.1, 2.0);
  snprintf(name, 200, "PtFakeP%s", aName);
  fPtFakeP = new TH1D(name, "fake protons Pt",                100, 0.1, 2.0);
  snprintf(name, 200, "PtRP%s", aName);
  fPtRP = new TH1D(name, "right pdg protons Pt",                 100, 0.1, 2.0);   

  snprintf(name, 200, "PtContP%s", aName);
  fPtContP = new TH1D(name, "contamination",                 100, 0.1, 2.0);   
  snprintf(name, 200, "PtContPi%s", aName);
  fPtContPi = new TH1D(name, "contamination",                 100, 0.1, 2.0);   
  snprintf(name, 200, "PtContMup%s", aName);
  fPtContMup = new TH1D(name, "contamination",                 100, 0.1, 2.0);   
  snprintf(name, 200, "PtContElp%s", aName);
  fPtContElp = new TH1D(name, "contamination",                 100, 0.1, 2.0);   
 
}

AliFemtoCutMonitorParticlePtPDG::AliFemtoCutMonitorParticlePtPDG(const AliFemtoCutMonitorParticlePtPDG &aCut):
  AliFemtoCutMonitor(),
  fPtPDG(0),ftpcHist(0),fPtGoodPi(0),fPtFakePi(0),fPtGoodK(0),fPtFakeK(0),
  fPtGoodP(0),fPtFakeP(0),fPtRPi(0),fPtRK(0),fPtRP(0),
  fPtContP(0),
  fPtContPi(0),
  fPtContMup(0),
  fPtContElp(0), 
  fMass(0.13957)
{
  // copy constructor
  if (fPtPDG) delete fPtPDG;
  fPtPDG = new TH2D(*aCut.fPtPDG);
  ftpcHist= new TH2D(*aCut.ftpcHist);
  fPtGoodPi= new TH1D(*aCut.fPtGoodPi);
  fPtFakePi= new TH1D(*aCut.fPtFakePi);
  fPtGoodK= new TH1D(*aCut.fPtGoodK);
  fPtFakeK= new TH1D(*aCut.fPtFakePi);
  fPtGoodP= new TH1D(*aCut.fPtGoodP);
  fPtFakeP= new TH1D(*aCut.fPtFakePi);
  fPtRPi= new TH1D(*aCut.fPtRPi);
  fPtRK= new TH1D(*aCut.fPtRK);
  fPtRP= new TH1D(*aCut.fPtRP);  
  
  fPtContP= new TH1D(*aCut.fPtContP);
  fPtContPi= new TH1D(*aCut.fPtContPi);
  fPtContMup= new TH1D(*aCut.fPtContMup);
  fPtContElp= new TH1D(*aCut.fPtContElp);
  
  fMass = aCut.fMass; 
}

AliFemtoCutMonitorParticlePtPDG::~AliFemtoCutMonitorParticlePtPDG()
{
  // Destructor
  delete fPtPDG;
  delete ftpcHist;
  delete fPtGoodPi;
  delete fPtFakePi;
  delete fPtGoodK;
  delete fPtFakeK;
  delete fPtGoodP;
  delete fPtFakeP; 
  delete fPtRPi;
  delete fPtRK;
  delete fPtRP;

  delete fPtContP;
  delete fPtContPi;
  delete fPtContMup;
  delete fPtContElp;

}

AliFemtoCutMonitorParticlePtPDG& AliFemtoCutMonitorParticlePtPDG::operator=(const AliFemtoCutMonitorParticlePtPDG& aCut)
{
  // assignment operator
  if (this == &aCut) 
    return *this;

  if (fPtPDG) delete fPtPDG;
  fPtPDG = new TH2D(*aCut.fPtPDG);
  
  if (ftpcHist) delete ftpcHist;
  ftpcHist = new TH2D(*aCut.ftpcHist);
  
   if (fPtGoodPi) delete fPtGoodPi;
  fPtGoodPi = new TH1D(*aCut.fPtGoodPi);
  
   if (fPtFakePi) delete fPtFakePi;
  fPtFakePi = new TH1D(*aCut.fPtFakePi);
  
   if (fPtRPi) delete fPtRPi;
  fPtRPi = new TH1D(*aCut.fPtRPi);
  
  if (fPtGoodK) delete fPtGoodK;
  fPtGoodK = new TH1D(*aCut.fPtGoodK);
  
   if (fPtFakeK) delete fPtFakeK;
  fPtFakeK = new TH1D(*aCut.fPtFakeK);
  
   if (fPtRK) delete fPtRK;
  fPtRK = new TH1D(*aCut.fPtRK);  
   
  if (fPtGoodP) delete fPtGoodP;
  fPtGoodP = new TH1D(*aCut.fPtGoodP);
  
   if (fPtFakeP) delete fPtFakeP;
  fPtFakeP = new TH1D(*aCut.fPtFakeP);
  
   if (fPtRP) delete fPtRP;
  fPtRP = new TH1D(*aCut.fPtRP);
 
   if (fPtContP) delete fPtContP;
  fPtContP = new TH1D(*aCut.fPtContP);

   if (fPtContPi) delete fPtContPi;
  fPtContPi = new TH1D(*aCut.fPtContPi);
  
     if (fPtContMup) delete fPtContMup;
  fPtContMup = new TH1D(*aCut.fPtContMup);

     if (fPtContElp) delete fPtContElp;
  fPtContElp = new TH1D(*aCut.fPtContElp);
  
  return *this;
}

AliFemtoString AliFemtoCutMonitorParticlePtPDG::Report(){ 
  // Prepare report from the execution
  string stemp = "*** AliFemtoCutMonitorParticlePtPDG report"; 
  AliFemtoString returnThis = stemp;
  return returnThis; 
}

void AliFemtoCutMonitorParticlePtPDG::Fill(const AliFemtoTrack* aTrack)
{
  // Fill in the monitor histograms with the values from the current track
  //  float tEnergy = ::sqrt(aTrack->P().mag2()+fMass*fMass);
  //  float tRapidity = 0.5*::log((tEnergy+aTrack->P().z())/(tEnergy-aTrack->P().z()));
  float tPt = ::sqrt((aTrack->P().x())*(aTrack->P().x())+(aTrack->P().y())*(aTrack->P().y()));
  //  float tEta = -TMath::Log(TMath::Tan(aTrack->P().theta()/2.0));
  //  float tPhi = aTrack->P().phi();
  float tP = ::sqrt((aTrack->P().z())*(aTrack->P().z())+(aTrack->P().x())*(aTrack->P().x())+(aTrack->P().y())*(aTrack->P().y()));;
  float dedx = aTrack->TPCsignalN();
  float w[10];
  w[0] = aTrack->PidProbElectron();
  w[1] = aTrack->PidProbMuon();
  w[2] = aTrack->PidProbPion();
  w[3] = aTrack->PidProbKaon(); 
  w[4] = aTrack->PidProbProton();
  
   Int_t pdg1=0;
  AliFemtoModelHiddenInfo *info = ( AliFemtoModelHiddenInfo *) aTrack->GetHiddenInfo();
  if(info)pdg1 = info->GetPDGPid();


//most probable particle  
  fPtGoodPi->Fill(tPt);
  fPtGoodP->Fill(tPt);
  fPtGoodK->Fill(tPt);

//contaminations 
  if (abs(pdg1)!=321)fPtFakeK->Fill(tPt);
  if (abs(pdg1)!=211)fPtFakePi->Fill(tPt);
  if (abs(pdg1)!=2212)fPtFakeP->Fill(tPt);

               
//contaminations for kaons 
  if (abs(pdg1)==2212)fPtContP->Fill(tPt);
  if (abs(pdg1)==211)fPtContPi->Fill(tPt);
  if (abs(pdg1)==13)fPtContMup->Fill(tPt);
  if (abs(pdg1)==11)fPtContElp->Fill(tPt);
             
  Float_t pdg=-1.0;
  if(abs(pdg1)==211){
   pdg=2.0;
   fPtRPi->Fill(tPt);
  }

  if(abs(pdg1)==321){
   pdg=3.0;
   fPtRK->Fill(tPt);
   }
  
  if(abs(pdg1)==2212){
   pdg=4.0;
   fPtRP->Fill(tPt);
   }
  
  if(abs(pdg1)==11)pdg=0.0; //+-electron
  if(abs(pdg1)==13)pdg=1.0;  //+-muon
  
  //cout<<"pdg from CutMonitor.."<<pdg1<<"pdg"<<pdg<<endl;
 
   
  fPtPDG->Fill(pdg, tPt);
  ftpcHist->Fill(tP,dedx);
  
 
}

void AliFemtoCutMonitorParticlePtPDG::Write()
{
  // Write out the relevant histograms
  
  fPtPDG->Write();
  ftpcHist->Write();
  fPtGoodPi->Write();
  fPtFakePi->Write();
  fPtGoodK->Write();
  fPtFakeK->Write();
  fPtGoodP->Write();
  fPtFakeP->Write(); 
  fPtRPi->Write();
  fPtRK->Write();
  fPtRP->Write();
  fPtContP->Write();
  fPtContPi->Write();
  fPtContMup->Write();
  fPtContElp->Write();
}

TList *AliFemtoCutMonitorParticlePtPDG::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fPtPDG);
  tOutputList->Add(ftpcHist);
  tOutputList->Add(fPtGoodPi);
  tOutputList->Add(fPtFakePi);
  tOutputList->Add(fPtGoodK);
  tOutputList->Add(fPtFakeK);
  tOutputList->Add(fPtGoodP);
  tOutputList->Add(fPtFakeP); 
  tOutputList->Add(fPtRPi);
  tOutputList->Add(fPtRK);
  tOutputList->Add(fPtRP);
 
  tOutputList->Add(fPtContP);
  tOutputList->Add(fPtContPi);
  tOutputList->Add(fPtContMup);
  tOutputList->Add(fPtContElp); 

  return tOutputList;
}
