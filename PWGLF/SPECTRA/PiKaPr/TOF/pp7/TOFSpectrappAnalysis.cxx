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

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                                                                       //
// Analysis for identified charged hadron spectra. TOF                   //
//                                                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TFile.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenCocktailEventHeader.h"

#include "AliPID.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliAnalysisFilter.h"
#include "AliCFContainer.h"
#include "AliCFFrame.h"
#include "AliCFGridSparse.h"
#include "TDatabasePDG.h"
#include "TROOT.h"
#include "TSystem.h"
#include "AliMultiplicity.h"

#include "AliLog.h"

#include "TOFSpectrappAnalysis.h"
//cambio
#include "AliESDUtils.h"
//fine cambio
ClassImp(TOFSpectrappAnalysis)

//________________________________________________________________________
TOFSpectrappAnalysis::TOFSpectrappAnalysis() 
: AliAnalysisTaskSE("TaskChargedHadron"), fESD(0), fESDtrackCuts(0), fESDpid(0x0),
  fMCtrue(0),
  fAlephParameters(),
  spdCorr(-1),
  multiplicity(0),
  ZPrimVertex(0),
  frun(0),
  frunOld(0),
  fLoadOCDB(kTRUE),
  correctTExp(kTRUE),
  calibrateESD(kTRUE), useT0TOF(kTRUE), timeResolution(96.), tuneTOFMC(kFALSE),
  fTreeTrack(0),
  //fTreeEv(0),
  hNumEv(0),
  tofCalib(0),
  t0maker(0),
  fDCAXY(-999), fDCAZ(-999), kselimppar(-999), fmatch(-999), fmom(-999), flength(-999), fsign(-999), ftimetof(-999), fexptimeel(-999), fexptimemu(-999), fexptimepi(-999), fexptimeka(-999), fexptimepr(-999), ftofchan(-999), feta(-999), fphi(-999), fmomtrasv(-999), t0track(-999), t0trackSigma(-999), sigmael(-999), sigmamu(-999), sigmapi(-999), sigmaka(-999), sigmapr(-999), TPCSignal(-999), TPCSigmaPI(-999), TPCSigmaKA(-999), TPCSigmaPR(-999), TOFlabel0(-999), TOFlabel1(-999), TOFlabel2(-999)    
{
  // default Constructor
  for(Int_t i=0; i<5;i++){r1[i]=-999;}
}


//________________________________________________________________________
TOFSpectrappAnalysis::TOFSpectrappAnalysis(const char *name) 
  : AliAnalysisTaskSE(name), fESD(0),fESDtrackCuts(0),fESDpid(0x0),
    fMCtrue(0),
    fAlephParameters(),
    spdCorr(-1),
    multiplicity(0),
    ZPrimVertex(-999),
    frun(0),
    frunOld(0),
    fLoadOCDB(kTRUE),
    correctTExp(kTRUE),
    calibrateESD(kTRUE), useT0TOF(kTRUE), timeResolution(96.), tuneTOFMC(kFALSE),
    fTreeTrack(0),
    //fTreeEv(0),
    hNumEv(0),
    tofCalib(0),
    t0maker(0),
    fDCAXY(-999), fDCAZ(-999), kselimppar(-999), fmatch(-999), fmom(-999), flength(-999), fsign(-999), ftimetof(-999), fexptimeel(-999), fexptimemu(-999), fexptimepi(-999), fexptimeka(-999), fexptimepr(-999), ftofchan(-999), feta(-999), fphi(-999), fmomtrasv(-999), t0track(-999), t0trackSigma(-999), sigmael(-999), sigmamu(-999), sigmapi(-999), sigmaka(-999), sigmapr(-999), TPCSignal(-999), TPCSigmaPI(-999), TPCSigmaKA(-999), TPCSigmaPR(-999), TOFlabel0(-999), TOFlabel1(-999), TOFlabel2(-999)    
{
  
  //
  // standard constructur which should be used
  //
  Printf("*** CONSTRUCTOR CALLED ****");
  // 
  /* real */
  for(Int_t i=0; i<5;i++){r1[i]=-999;}


  fAlephParameters[0] = 0.0283086;
  fAlephParameters[1] = 2.63394e+01;
  fAlephParameters[2] = 5.04114e-11;
  fAlephParameters[3] = 2.12543e+00;
  fAlephParameters[4] = 4.88663e+00;
  //
  // initialize PID object
  //

  
  //tofCalib = new AliTOFcalib();

  //cambio

  //fESDpid = new AliESDpid();

  //fine cambio
  //t0maker = new AliTOFT0maker(fESDpid, tofCalib); 
  //
  // create track cuts
  //
  AliESDtrackCuts* ESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts");
  fESDtrackCuts = (AliESDtrackCuts*) ESDtrackCuts->GetStandardITSTPCTrackCuts2010(kTRUE);

  //fESDpid->GetTPCResponse().SetBetheBlochParameters(fAlephParameters[0],fAlephParameters[1],fAlephParameters[2],fAlephParameters[3],fAlephParameters[4]);

  //
  //Initialize();
  // Output slot #0 writes into a TList container
  //DefineOutput(1, TTree::Class());
  DefineOutput(1, TTree::Class());
  DefineOutput(2, TH1D::Class());

}


//________________________________________________________________________
Int_t TOFSpectrappAnalysis::Mult()
{
  AliESDtrackCuts* esdTrackCutsMult = new AliESDtrackCuts;
  
  // TPC
  esdTrackCutsMult->SetMinNClustersTPC(70);//ok
  esdTrackCutsMult->SetMaxChi2PerClusterTPC(4);//ok
  esdTrackCutsMult->SetAcceptKinkDaughters(kFALSE);//ok
  esdTrackCutsMult->SetRequireTPCRefit(kTRUE);//ok
  // ITS
  esdTrackCutsMult->SetRequireITSRefit(kTRUE);//ok
  esdTrackCutsMult->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					     AliESDtrackCuts::kAny);//ok
  
  esdTrackCutsMult->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
  
  esdTrackCutsMult->SetMaxDCAToVertexZ(2);//ok
  
  esdTrackCutsMult->SetDCAToVertex2D(kFALSE);//ok
  esdTrackCutsMult->SetRequireSigmaToVertex(kFALSE);//ok
  
  esdTrackCutsMult->SetEtaRange(-0.8,+0.8);
  esdTrackCutsMult->SetPtRange(0.15, 1e10);
  
  // Here is what we finally want:
  Int_t multipli = esdTrackCutsMult->CountAcceptedTracks(fESD);
  delete esdTrackCutsMult;
  return multipli;
  
 }


//________________________________________________________________________
void TOFSpectrappAnalysis::UserCreateOutputObjects() 
{
  // Create histograms
  // Called once
  
  OpenFile(1);

  fTreeTrack = new TTree("TreeTrack","Track Properties");
  fTreeTrack->Branch("spdCorr",&spdCorr,"spdCorr/F");  
  fTreeTrack->Branch("multiplicity",&multiplicity,"multiplicity/I");
  fTreeTrack->Branch("DCAXY",&fDCAXY,"fDCAXY/F");  
  fTreeTrack->Branch("DCAZ",&fDCAZ,"fDCAZ/F");  
  //fTreeTrack->Branch("ZPrimVertex",&ZPrimVertex,"ZPrimVertex/D");  
  //fTreeTrack->Branch("kselimppar",&kselimppar,"kselimppar/I");  
  fTreeTrack->Branch("fmatch",&fmatch,"fmatch/I");
  fTreeTrack->Branch("fmom",&fmom,"fmom/D");
  fTreeTrack->Branch("length",&flength,"length/D");
  fTreeTrack->Branch("sign",&fsign,"sign/I");
  fTreeTrack->Branch("timetof",&ftimetof,"timetof/D");
  fTreeTrack->Branch("exptimeel",&fexptimeel,"exptimeel/D");
  fTreeTrack->Branch("exptimemu",&fexptimemu,"exptimemu/D");
  fTreeTrack->Branch("exptimepi",&fexptimepi,"exptimepi/D");
  fTreeTrack->Branch("exptimeka",&fexptimeka,"exptimeka/D");
  fTreeTrack->Branch("exptimepr",&fexptimepr,"exptimepr/D");
  fTreeTrack->Branch("tofchan",&ftofchan,"tofchan/I"); 
  fTreeTrack->Branch("eta",&feta,"eta/D");
  fTreeTrack->Branch("phi",&fphi,"phi/D");
  fTreeTrack->Branch("momtrasv",&fmomtrasv,"momtrasv/D");
  fTreeTrack->Branch("t0trackSigma",&t0trackSigma,"t0trackSigma/F");
  fTreeTrack->Branch("t0track",&t0track,"t0track/F");
  fTreeTrack->Branch("sigmael",&sigmael,"sigmael/D");
  fTreeTrack->Branch("sigmamu",&sigmamu,"sigmamu/D");
  fTreeTrack->Branch("sigmapi",&sigmapi,"sigmapi/D");
  fTreeTrack->Branch("sigmaka",&sigmaka,"sigmaka/D");
  fTreeTrack->Branch("sigmapr",&sigmapr,"sigmapr/D");
 //  //fTreeTrack->Branch("TPCsignal",&TPCSignal,"TPCsignal/D");
//   //fTreeTrack->Branch("TPCSigmaPI",&TPCSigmaPI,"TPCSigmaPI/F");
//   //fTreeTrack->Branch("TPCSigmaKA",&TPCSigmaKA,"TPCSigmaKA/F");
//   //fTreeTrack->Branch("TPCSigmaPR",&TPCSigmaPR,"TPCSigmaPR/F");
//   fTreeTrack->Branch("r10",&r1[0],"r10/D");
//   fTreeTrack->Branch("r11",&r1[1],"r11/D");
//   fTreeTrack->Branch("r12",&r1[2],"r12/D");
//   fTreeTrack->Branch("r13",&r1[3],"r13/D");
//   fTreeTrack->Branch("r14",&r1[4],"r14/D");
//   // fTreeTrack->Branch("TOFlabel0",&TOFlabel0,"TOFlabel0/I");
// //   fTreeTrack->Branch("TOFlabel1",&TOFlabel1,"TOFlabel1/I");
// //   fTreeTrack->Branch("TOFlabel2",&TOFlabel2,"TOFlabel2/I");

 
  PostData(1, fTreeTrack );

  // fTreeEv = new TTree("TreeEv","Event Properties");
//   fTreeEv->Branch("multiplicity",&multiplicity,"multiplicity/I");  
//   fTreeEv->Branch("spdCorr",&spdCorr,"spdCorr/F");  
//   //fTreeEv->Branch("ZPrimVertex",&ZPrimVertex,"ZPrimVertex/D"); 
  

  OpenFile(2);
  hNumEv=new TH1D("NemEv","NemEv",4,1,5);
  PostData(2, hNumEv);

}

//________________________________________________________________________
void TOFSpectrappAnalysis::UserExec(Option_t *) 
{
  //
  // main event loop
  //
 
  
  multiplicity=-999, ZPrimVertex=-999, frun=-999, spdCorr=-1.0;

  fESD = dynamic_cast<AliESDEvent*>( InputEvent() );
  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }
  
  if (!fESDtrackCuts) {
    Printf("ERROR: fESDtrackCuts not available");
    return;
  }

  hNumEv->Fill(1);

  //cambio
  AliESDInputHandler* esdH = dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (esdH)
    fESDpid = esdH->GetESDpid();
  //fine cambio
  //
  // check if event is selected by physics selection class
  //
  Bool_t isSelected = kFALSE;
  isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()& AliVEvent::kMB);

  if (!isSelected) {
    return;
  }
  
  hNumEv->Fill(2);

  //
  // monitor vertex position
  //
  const AliESDVertex *vertex = fESD->GetPrimaryVertexTracks(); //! Primary vertex estimated using ESD tracks
  if(vertex->GetNContributors()<1) { // # of tracklets/tracks used for the estimate
    // SPD vertex
    vertex = fESD->GetPrimaryVertexSPD(); //! Primary vertex estimated by the SPD
    if(vertex->GetNContributors()<1) vertex = 0x0;
  }  

  if (!vertex) {
    return;
  } 
  
  hNumEv->Fill(3);

  if (vertex) {ZPrimVertex=vertex->GetZv();}
  if(TMath::Abs(ZPrimVertex)>10){return;}
  
  hNumEv->Fill(4);
  
  multiplicity=Mult();
  
  //multiplicity as defined by Marek
  const AliMultiplicity *mult = fESD->GetMultiplicity();
  Float_t nClusters[6]={0.0,0.0,0.0,0.0,0.0,0.0};
  for(Int_t ilay=0; ilay<6; ilay++)
    {
      nClusters[ilay] = (Float_t)mult->GetNumberOfITSClusters(ilay);
    }
  //cambio
  spdCorr = AliESDUtils::GetCorrSPD2(nClusters[1],vertex->GetZ());
  //fine cambio
  //fTreeEv->Fill();
  

  //TOF settings done in the TOF tender
  
  frun = fESD->GetRunNumber();
 //  if(frun==frunOld){fLoadOCDB=kFALSE;}else {fLoadOCDB=kTRUE;}
//   //if (tuneTOFMC) calibrateESD = kFALSE;
//   Double_t *T0TOF;
//   if(fLoadOCDB){
//   AliCDBManager *cdb = AliCDBManager::Instance();
//   cdb->SetDefaultStorage("alien://folder=/alice/data/2010/OCDB");
//   //cdb->SetDefaultStorage("raw://");
//   cdb->SetRun(frun);}
  
//   /* init TOF calibration */
//   if (correctTExp)
//     tofCalib->SetCorrectTExp(kTRUE);
//   tofCalib->Init();
  
//   /* init TOF T0-maker */
//   t0maker->SetTimeResolution(timeResolution);
  
//   /* calibrate ESD */
//   if (calibrateESD)
//     tofCalib->CalibrateESD(fESD);
   
//   T0TOF=t0maker->ComputeT0TOF(fESD);// calcola il t0 solo col tof in 10 bin di pt e lo setta in AliTOFPidResponse 

//   //scrive i valori precedenti nel TOFHeader
//   t0maker->WriteInESD(fESD); 

//   //setta T0_TOF in AliTOFPidResponse ovvero i valori settati in AliTOFHeader. Se non ci sono setta T0spread
//   fESDpid->SetTOFResponse(fESD,AliESDpid::kTOF_T0); 


//   fESDpid->MakePID(fESD,kFALSE); //calcola la sigma e le gi in più definisce la flag kTOFmismatch
  

  //
  // track loop
  //
  
  //const Float_t kNsigmaCut = 3;
  //const Float_t k2sigmaCorr = 1/(0.5*(TMath::Erf(kNsigmaCut/sqrt(2))-TMath::Erf(-kNsigmaCut/sqrt(2))))/*1/0.9545*/;
  //Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z for the vertex cut
  
  for (Int_t i=0;i<fESD->GetNumberOfTracks();++i) {
    
    AliESDtrack *track =fESD->GetTrack(i); 
    if (!track){continue;}
    // start TOF analysis
    
    fDCAXY=-999, fDCAZ=-999, kselimppar=-999, fmatch=-999, fmom=-999, flength=-999, fsign=-999, ftimetof=-999, fexptimeel=-999, fexptimeel=-999, fexptimepi=-999, fexptimeka=-999, fexptimepr=-999, ftofchan=-999, feta=-999, fphi=-999, fmomtrasv=-999, sigmael=-999, sigmamu=-999, sigmapi=-999, t0track=-999, t0trackSigma=-999, sigmaka=-999, sigmapr=-999, TPCSignal=-999, TPCSigmaPI=-999, TPCSigmaKA=-999, TPCSigmaPR=-999, r1[0]=-999,r1[1]=-999,r1[2]=-999,r1[3]=-999,r1[4]=-999;
    
    if (!fESDtrackCuts->AcceptTrack(track)) {continue;}
    if (SelectOnImpPar(track)) {kselimppar=1;}else {kselimppar=0;}
    if ((track->GetStatus()&AliESDtrack::kTOFout)==0) {continue;}
    if ((track->GetStatus()&AliESDtrack::kTIME)==0) {continue;}  
    if (!(track->GetStatus()&AliESDtrack::kTOFmismatch)==0) {fmatch=0;}else {fmatch=1;}  
    track->GetImpactParameters(fDCAXY, fDCAZ);
    track->GetTOFpid(r1);
    fmom=track->GetP();  // This function returns the track momentum
    flength=track->GetIntegratedLength();
    fsign=track->GetSign();
    ftimetof=track->GetTOFsignal();
    Double_t inttime[5]; 
    track->GetIntegratedTimes(inttime);// Returns the array with integrated times for each particle hypothesis
    fexptimeel=inttime[0];
    fexptimemu=inttime[1];
    fexptimepi=inttime[2];
    fexptimeka=inttime[3];
    fexptimepr=inttime[4];  
    ftofchan=track->GetTOFCalChannel(); // Channel Index of the TOF Signal
    feta=track->Eta();// return pseudorapidity return -TMath::Log(TMath::Tan(0.5 * Theta())); 
    fphi=track->Phi()* 180 / TMath::Pi();// Returns the azimuthal angle of momentum  0 <= phi < 2*pi
    fmomtrasv=track->Pt(); 
    t0track = fESDpid->GetTOFResponse().GetStartTime(fmom); // T0best time
    t0trackSigma = fESDpid->GetTOFResponse().GetStartTimeRes(fmom); // T0best time
    Double_t sigma[5];
    for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++){
      sigma[ipart] = fESDpid->GetTOFResponse().GetExpectedSigma(fmom, inttime[ipart], AliPID::ParticleMass(ipart));}
    sigmael=sigma[0];
    sigmamu=sigma[1];
    sigmapi=sigma[2];
    sigmaka=sigma[3];
    sigmapr=sigma[4];
    Int_t TOFLAB[3];
    track->GetTOFLabel(TOFLAB);
    TOFlabel0=TOFLAB[0];
    TOFlabel1=TOFLAB[1];
    TOFlabel2=TOFLAB[2];
    //TPCSignal = track->GetTPCsignal();
    //TPCSigmaPI=fESDpid->NumberOfSigmasTPC(track,AliPID::kPion);
    //TPCSigmaKA=fESDpid->NumberOfSigmasTPC(track,AliPID::kKaon);
    //TPCSigmaPR=fESDpid->NumberOfSigmasTPC(track,AliPID::kProton);
    fTreeTrack->Fill(); 
    
    /*
      kTOFout = TOF matching
      kTIME = good integrated time
      100000 > track->GetTOFsignal() > 12000 = TOF time reasanble range
      tracklength > 365 = should be greater than the TOF radius (370 cm)
    */
    
    
    
  
    
    
  } // end of track loop
  
  frunOld=frun;
  
  
  // Post output data  
  PostData(1, fTreeTrack);
  //PostData(1, fTreeEv);
  PostData(2, hNumEv);

}      


//________________________________________________________________________
void TOFSpectrappAnalysis::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  Printf("*** CONSTRUCTOR CALLED ****");
  fTreeTrack = dynamic_cast<TTree*> (GetOutputData(1));
  if (!fTreeTrack) {
    Printf("ERROR: fTreeTrack not available");
    return;   
  } 
  
  system("touch ok.job");
}


//________________________________________________________________________
Bool_t TOFSpectrappAnalysis::SelectOnImpPar(AliESDtrack* t) {
  //
  // cut on transverse impact parameter
  //
  Float_t d0z0[2],covd0z0[3];
  t->GetImpactParameters(d0z0,covd0z0);
  Float_t sigma= 0.0050+0.0060/TMath::Power(t->Pt(),0.9);
  Float_t d0max = 7.*sigma;
  //
  Float_t sigmaZ = 0.0146+0.0070/TMath::Power(t->Pt(),1.114758);
  if (t->Pt() > 1) sigmaZ = 0.0216;
  Float_t d0maxZ = 5.*sigmaZ;
  //
  if(TMath::Abs(d0z0[0]) < d0max && d0z0[1] < d0maxZ) return kTRUE;
  return kFALSE;
}