////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCutMonitorParticlePtPDGV0 - the cut monitor for particles to study //
// the difference between reconstructed V0s and MC V0s                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitorParticlePtPDGV0.h"
#include "AliFemtoModelHiddenInfo.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TMath.h>

AliFemtoCutMonitorParticlePtPDGV0::AliFemtoCutMonitorParticlePtPDGV0():
  fPtPDG(0),ftpcHist(0),fPtMostProbable(0),
  fPtFakeLambdas(0), fFakeProtonDaughters(0), fFakeAntiProtonDaughters(0), fFakePionPlusDaughters(0), fFakePionMinusDaughters(0),
  //fPtRPi(0),fPtRK(0),fPtRP(0),
  fPtV0(0), fPtPosProton(0),
  fPtNegProton(0),
  fPtPosPion(0),
  fPtNegPion(0),
  fMass(0.13957)
{
  // Default constructor
}

AliFemtoCutMonitorParticlePtPDGV0::AliFemtoCutMonitorParticlePtPDGV0(const char *aName, float aMass):
  AliFemtoCutMonitor(),
  fPtPDG(0),ftpcHist(0),fPtMostProbable(0),
  fPtFakeLambdas(0), fFakeProtonDaughters(0), fFakeAntiProtonDaughters(0), fFakePionPlusDaughters(0), fFakePionMinusDaughters(0),
  //fPtRPi(0),fPtRK(0),fPtRP(0),
  fPtV0(0), fPtPosProton(0),
  fPtNegProton(0),
  fPtPosPion(0),
  fPtNegPion(0),
  fMass(aMass)
{
  // Normal constructor
  char name[200];
  snprintf(name, 200, "PtPDG%s", aName);
  fPtPDG = new TH2D(name, "PDG vs Pt", 10, 0.0, 5.0, 8000, -4000, 4000);
  snprintf(name, 200, "tpcHist%s", aName);
  ftpcHist=new TH2D(name,"TPC dE/dX vs momentum",100,0.1,2.7,100,0.,6.);
  snprintf(name, 200, "PtOfV0s%s", aName);
  fPtMostProbable = new TH1D(name, "Pt of V0s",      250, 0.1, 5.0);
  snprintf(name, 200, "PtFakeLambdas%s", aName);
  fPtFakeLambdas = new TH1D(name, "fake lambdas Pt",              250, 0.1, 5.0);
  snprintf(name, 200, "PtFakePionPlus%s", aName);
  fFakePionPlusDaughters = new TH1D(name, "fake pion+ daughters",              250, 0.1, 5.0);
  snprintf(name, 200, "PtFakeProtonPlus%s", aName);
  fFakeProtonDaughters = new TH1D(name, "fake proton daughters",              250, 0.1, 5.0);
  snprintf(name, 200, "PtFakePionMinus%s", aName);
  fFakePionMinusDaughters = new TH1D(name, "fake pion- daughters",              250, 0.1, 5.0);
  snprintf(name, 200, "PtFakeProtonPlus%s", aName);
  fFakeAntiProtonDaughters = new TH1D(name, "fake anti-proton daughters",              250, 0.1, 5.0);



  // fPtRPi = new TH1D(name, "right pdg pions Pt",               250, 0.1, 5.0);
  // snprintf(name, 200, "PtFakeK%s", aName);
  // snprintf(name, 200, "PtRK%s", aName);
  // fPtRK = new TH1D(name, "right pdg kaons Pt",
  // snprintf(name, 200, "PtRP%s", aName);
  // fPtRP = new TH1D(name, "right pdg protons Pt",                 250, 0.1, 5.0);   

  snprintf(name, 200, "PtRealLambdas%s", aName);
  fPtV0 = new TH1D(name, "Pt of real Lambdas in the sample", 250, 0.1, 5.0);   
  snprintf(name, 200, "PtRealPosP%s", aName);
  fPtPosProton = new TH1D(name, "Pt of real protons in the sample",  250, 0.1, 5.0);  
  snprintf(name, 200, "PtRealNegP%s", aName);
  fPtNegProton = new TH1D(name, "Pt of real antiprotons in the sample", 250, 0.1, 5.0);   
  snprintf(name, 200, "PtRealPosPi%s", aName);
  fPtPosPion = new TH1D(name, "Pt of real pion+ in the sample", 250, 0.1, 5.0);   
  snprintf(name, 200, "PtRealNegPi%s", aName);
  fPtNegPion = new TH1D(name, "Pt of real pion- in the sample", 250, 0.1, 5.0);   
 
}

AliFemtoCutMonitorParticlePtPDGV0::AliFemtoCutMonitorParticlePtPDGV0(const AliFemtoCutMonitorParticlePtPDGV0 &aCut):
  AliFemtoCutMonitor(),
  fPtPDG(0),ftpcHist(0),fPtMostProbable(0),
  fPtFakeLambdas(0), fFakeProtonDaughters(0), fFakeAntiProtonDaughters(0), fFakePionPlusDaughters(0), fFakePionMinusDaughters(0),
  //fPtRPi(0),fPtRK(0),fPtRP(0),
  fPtV0(0), fPtPosProton(0),
  fPtNegProton(0),
  fPtPosPion(0),
  fPtNegPion(0), 
  fMass(0.13957)
{
  // copy constructor
  if (fPtPDG) delete fPtPDG;
  fPtPDG = new TH2D(*aCut.fPtPDG);
  ftpcHist= new TH2D(*aCut.ftpcHist);
  fPtMostProbable= new TH1D(*aCut.fPtMostProbable);
  fPtFakeLambdas= new TH1D(*aCut.fPtFakeLambdas);
  fFakeProtonDaughters= new TH1D(*aCut.fFakeProtonDaughters);
  fFakeAntiProtonDaughters= new TH1D(*aCut.fFakeAntiProtonDaughters);
  fFakePionPlusDaughters= new TH1D(*aCut.fFakePionPlusDaughters);
  fFakePionMinusDaughters= new TH1D(*aCut.fFakePionMinusDaughters);

  // fPtRPi= new TH1D(*aCut.fPtRPi);
  // fPtRK= new TH1D(*aCut.fPtRK);
  // fPtRP= new TH1D(*aCut.fPtRP);  




  fPtV0= new TH1D(*aCut.fPtV0);
  fPtPosProton= new TH1D(*aCut.fPtPosProton);
  fPtNegProton= new TH1D(*aCut.fPtNegProton);
  fPtPosPion= new TH1D(*aCut.fPtPosPion);
  fPtNegPion= new TH1D(*aCut.fPtNegPion);
  
  fMass = aCut.fMass; 
}

AliFemtoCutMonitorParticlePtPDGV0::~AliFemtoCutMonitorParticlePtPDGV0()
{
  // Destructor
  delete fPtPDG;
  delete ftpcHist;
  delete fPtMostProbable;
  // delete fPtRPi;
  // delete fPtRK;
  // delete fPtRP;

  delete fPtFakeLambdas;
  delete fFakeProtonDaughters;
  delete fFakePionPlusDaughters; 
  delete fFakeAntiProtonDaughters;
  delete fFakePionMinusDaughters;

  delete fPtV0;
  delete fPtPosProton;
  delete fPtNegProton;
  delete fPtPosPion;
  delete fPtNegPion;

}

AliFemtoCutMonitorParticlePtPDGV0& AliFemtoCutMonitorParticlePtPDGV0::operator=(const AliFemtoCutMonitorParticlePtPDGV0& aCut)
{
  // assignment operator
  if (this == &aCut) 
    return *this;

  if (fPtPDG) delete fPtPDG;
  fPtPDG = new TH2D(*aCut.fPtPDG);
  
  if (ftpcHist) delete ftpcHist;
  ftpcHist = new TH2D(*aCut.ftpcHist);
  
   if (fPtFakeLambdas) delete fPtFakeLambdas;
  fPtFakeLambdas = new TH1D(*aCut.fPtFakeLambdas);
  
   if (fFakeProtonDaughters) delete fFakeProtonDaughters;
   fFakeProtonDaughters = new TH1D(*aCut.fFakeProtonDaughters);
    
   if (fFakeAntiProtonDaughters) delete fFakeAntiProtonDaughters;
   fFakeAntiProtonDaughters = new TH1D(*aCut.fFakeAntiProtonDaughters);
     
   if (fFakePionPlusDaughters) delete fFakePionPlusDaughters;
   fFakePionPlusDaughters = new TH1D(*aCut.fFakePionPlusDaughters);
      
   if (fFakePionMinusDaughters) delete fFakePionMinusDaughters;
   fFakePionMinusDaughters = new TH1D(*aCut.fFakePionMinusDaughters); 

  //  if (fPtRPi) delete fPtRPi;
  // fPtRPi = new TH1D(*aCut.fPtRPi);
  
  //  if (fPtRK) delete fPtRK;
  // fPtRK = new TH1D(*aCut.fPtRK);  

  //  if (fPtRP) delete fPtRP;
  // fPtRP = new TH1D(*aCut.fPtRP);
 
   if (fPtV0) delete fPtV0;
  fPtV0 = new TH1D(*aCut.fPtV0);
 
   if (fPtPosProton) delete fPtPosProton;
  fPtPosProton = new TH1D(*aCut.fPtPosProton);

   if (fPtNegProton) delete fPtNegProton;
  fPtNegProton = new TH1D(*aCut.fPtNegProton);
  
     if (fPtPosPion) delete fPtPosPion;
  fPtPosPion = new TH1D(*aCut.fPtPosPion);

     if (fPtNegPion) delete fPtNegPion;
  fPtNegPion = new TH1D(*aCut.fPtNegPion);
  
  return *this;
}

AliFemtoString AliFemtoCutMonitorParticlePtPDGV0::Report(){ 
  // Prepare report from the execution
  string stemp = "*** AliFemtoCutMonitorParticlePtPDGV0 report"; 
  AliFemtoString returnThis = stemp;
  return returnThis; 
}

void AliFemtoCutMonitorParticlePtPDGV0::Fill(const AliFemtoV0* aV0)
{
  
  // Fill in the monitor histograms with the values from the current track
  //  float tEnergy = ::sqrt(aV0->P().mag2()+fMass*fMass);
  //  float tRapidity = 0.5*::log((tEnergy+aV0->P().z())/(tEnergy-aV0->P().z()));
  float tPt = aV0->PtV0();//::sqrt((aV0->P().x())*(aV0->P().x())+(aV0->P().y())*(aV0->P().y()));
  float tPtPos = aV0->PtPos();
  float tPtNeg = aV0->PtNeg();
  //  float tEta = -TMath::Log(TMath::Tan(aV0->P().theta()/2.0));
  //  float tPhi = aV0->P().phi();
  //float tP = ::sqrt((aV0->P().z())*(aV0->P().z())+(aV0->P().x())*(aV0->P().x())+(aV0->P().y())*(aV0->P().y()));;
  // float dedx = aV0->TPCsignalN();
  // float w[10];
  // w[0] = aV0->PidProbElectron();
  // w[1] = aV0->PidProbMuon();
  // w[2] = aV0->PidProbPion();
  // w[3] = aV0->PidProbKaon(); 
  // w[4] = aV0->PidProbProton();
  
  Int_t pdg1=0, pdgPosDaughter=0, pdgNegDaughter=0;
  AliFemtoModelHiddenInfo *info = ( AliFemtoModelHiddenInfo *) aV0->GetHiddenInfo();
  if(info)pdg1 = info->GetPDGPid();
  if(info)pdgPosDaughter=info->GetPDGPidPos();
  if(info)pdgNegDaughter=info->GetPDGPidNeg();
  

//most probable particle  
  fPtMostProbable->Fill(tPt);

//contaminations 
  if (abs(pdg1)!=3122)fPtFakeLambdas->Fill(tPt);
  if (pdgPosDaughter!=2212)fFakeProtonDaughters->Fill(tPtPos);
  if (pdgPosDaughter!=211)fFakePionPlusDaughters->Fill(tPtPos);
  if (pdgNegDaughter!=-2212)fFakeAntiProtonDaughters->Fill(tPtNeg);
  if (pdgNegDaughter!=-211)fFakePionMinusDaughters->Fill(tPtNeg);

               
//spectra of different particle species we've got
  if (abs(pdg1)==3122)fPtV0->Fill(tPt);
  if (abs(pdgNegDaughter)==2212)fPtNegProton->Fill(tPtNeg);
  if (abs(pdgPosDaughter)==2212)fPtPosProton->Fill(tPtPos);
  if (abs(pdgPosDaughter)==211)fPtPosPion->Fill(tPtPos);
  if (abs(pdgNegDaughter)==211)fPtNegPion->Fill(tPtNeg);

  Float_t pdg=pdg1;
  //if(abs(pdg1)==3122)pdg=0.0; //lambda
  //if(abs(pdg1)==211) pdg=1.0; //pion
  //if(abs(pdg1)==2212)pdg=2.0; //proton
  
  //cout<<"pdg from CutMonitor.."<<pdg1<<"pdg"<<pdg<<endl;
 
   
  fPtPDG->Fill(pdg, tPt);
  // ftpcHist->Fill(tP,dedx);
  
 
}

void AliFemtoCutMonitorParticlePtPDGV0::Write()
{
  // Write out the relevant histograms
  
  fPtPDG->Write();
  //ftpcHist->Write();
  fPtMostProbable->Write();
  fPtFakeLambdas->Write();
  fFakeProtonDaughters->Write();
  fFakePionPlusDaughters->Write();
  fFakeAntiProtonDaughters->Write();
  fFakePionMinusDaughters->Write();
  // fPtRPi->Write();
  // fPtRK->Write();
  // fPtRP->Write();
  fPtV0->Write();
  fPtPosProton->Write();
  fPtNegProton->Write();
  fPtPosPion->Write();
  fPtNegPion->Write();
}

TList *AliFemtoCutMonitorParticlePtPDGV0::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fPtPDG);
  //tOutputList->Add(ftpcHist);
  tOutputList->Add(fPtMostProbable);
  tOutputList->Add(fPtFakeLambdas);
  tOutputList->Add(fFakeProtonDaughters);
  tOutputList->Add(fFakePionPlusDaughters);
  tOutputList->Add(fFakeAntiProtonDaughters);
  tOutputList->Add(fFakePionMinusDaughters);
 
  tOutputList->Add(fPtV0);
  tOutputList->Add(fPtPosProton);
  tOutputList->Add(fPtNegProton);
  tOutputList->Add(fPtPosPion);
  tOutputList->Add(fPtNegPion); 

  return tOutputList;
}
