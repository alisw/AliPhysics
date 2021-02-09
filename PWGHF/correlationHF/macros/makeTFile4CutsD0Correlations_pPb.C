#include <Riostream.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TParameter.h>

//Use:
//Set hard coded commentet with //set this!!
// root[] .L makeInputD0tasks.C++
// root[] makeInputAliAnalysisTaskSED0Correlations()
//similar macros for the other D mesons

//Author: Fabio Colamaria, fabio.colamaria@ba.infn.it

//macro to make a .root file which contains an AliRDHFCutsD0toKpi for AliAnalysisTaskSED0Mass task

/***************************************************************************/
/***************************************************************************/
//  
//  Selections for pPb 2016 cent-integrated preliminaries (SQM2017)
//
//  *** D0/Event CUTS ***
//  - Optimized cuts (std2013 pPb + tighter cosThPoint + NormLxy + topomatic)
//  - Pileup rejection with 5 contributors, <0.8 cm (F.Prino)
//  - ZNA centrality estimator (though 0-100 centrality used)
//  - Reference pT bins (0-0.5-1-2-3-4-5-6-7-8-10-12-16-24-36)
//  
//  ----> For MC: no pileup selection; don't use centrality
//
//  *** ASSOC. TRACK CUTS ***
//  - Filterbit0
//  - 2 ITS cls, no ITSrefit
//  - TPC cuts using crossed rows
//  - DCAxy = DCAz = 1 cm (not 0.25cm)
//  - SB ranges in reference pT bins
//  - Optimized pools: -10,-1.5,3.5,10 cm; 0,35,55,500 trkl
//  
/***************************************************************************/
/***************************************************************************/

void makeInputAliAnalysisTaskSED0Correlations_pPb(Bool_t isOnData = kTRUE){

//____________________________________________________

  //Set Centrality
  Float_t minc=0,maxc=100;

// Cuts for D0 cuts

  AliRDHFCutsD0toKpi* RDHFD0Corr=new AliRDHFCutsD0toKpi();
  RDHFD0Corr->SetName("D0toKpiCuts");
  RDHFD0Corr->SetTitle("Cuts for D0 analysis");

  // PILE UP REJECTION
  if(isOnData) {
    RDHFD0Corr->SetOptPileup(AliRDHFCuts::kRejectPileupEvent);  	   //per DATI (spegni per MC)		
    RDHFD0Corr->ConfigurePileupCuts(5,0.8);  //-> Overloaded by mult dep pileup (F.Prino) //per DATI (spegni per MC) 
    //RDHFD0Corr->SetUseMultDepPileupCut(kTRUE);  	   //NON USATO PIU! //per DATI (spegni per MC)		 
  }

  //Event cuts
  RDHFD0Corr->SetMinVtxContr(1);
  RDHFD0Corr->SetMaxVtxZ(10.);
  RDHFD0Corr->SetCutOnzVertexSPD(2); 

  RDHFD0Corr->SetSelectCandTrackSPDFirst(kTRUE, 4); //Cristina 2016...

  //Trigger selection
  RDHFD0Corr->SetTriggerClass("");
  if(isOnData) RDHFD0Corr->SetTriggerMask(AliVEvent::kINT7);
  else RDHFD0Corr->SetTriggerMask(AliVEvent::kINT7); //A&E said to prefer kINT7 to kMB also in pPb MC

  //Quality tracks for daughters
  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  //esdTrackCuts->SetMinNClustersITS(4); // default is 5
  //esdTrackCuts->SetMinNClustersTPC(120);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny); // default is kBoth, otherwise kAny
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetPtRange(0.3,1.e10);

  RDHFD0Corr->AddTrackCuts(esdTrackCuts);

  //D0 selection topological cuts
  const Int_t nptbins =14;
  const Double_t ptmax = 9999.;
  const Int_t nvars=11;
  Float_t ptbins[nptbins+1];
  ptbins[0]=0.;
  ptbins[1]=0.5;	
  ptbins[2]=1.;
  ptbins[3]=2.;
  ptbins[4]=3.;
  ptbins[5]=4.;
  ptbins[6]=5.;
  ptbins[7]=6.;
  ptbins[8]=7.;
  ptbins[9]=8.;
  ptbins[10]=10.;
  ptbins[11]=12.;
  ptbins[12]=16.;
  ptbins[13]=24.;
  ptbins[14]=ptmax;

  RDHFD0Corr->SetGlobalIndex(nvars,nptbins);
  RDHFD0Corr->SetPtBins(nptbins+1,ptbins);
												 //        dca    ThSt ptk ptp     d0k       d0p         d02      thp
  Float_t cutsMatrixD0toKpiStand[nptbins][nvars]={{0.400,300.*1E-4,0.8,0.5,0.5,1000.*1E-4,1000.*1E-4,-50000.*1E-8,0.95,0.,5.},   /*BIN 0*/ /* pt<0.5*/ 
                                                  {0.400,300.*1E-4,0.8,0.5,0.5,1000.*1E-4,1000.*1E-4,-50000.*1E-8,0.95,0.,5.},   /*BIN 1*/ /* 0.5<pt<1*/
                                                  {0.400,300.*1E-4,0.8,0.5,0.5,1000.*1E-4,1000.*1E-4,-35000.*1E-8,0.95,0.,5.}, /*BIN 2*/ /* 1<pt<2 */
                                                  {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-30000.*1E-8,0.95,0.,5.}, /*BIN 3*/ /* 2<pt<3 */
                                                  {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-30000.*1E-8,0.95,0.,5.}, /*BIN 4*/ /* 3<pt<4 */
                                                  {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-15000.*1E-8,0.95,0.,5.}, /*BIN 5*/ /* 4<pt<5 */
                                                  {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-10000.*1E-8,0.95,0.,4.}, /*BIN 6*/ /* 5<pt<6 */
                                                  {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-8000.*1E-8,0.95,0.,4.},  /*BIN 7*/ /* 6<pt<7 */
                                                  {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-8000.*1E-8,0.95,0.,4.},  /*BIN 8*/ /* 7<pt<8 */
                                                  {0.400,300.*1E-4,0.9,0.7,0.7,1000.*1E-4,1000.*1E-4,-5000.*1E-8,0.95,0.,3.},  /*BIN 9*/ /* 8<pt<10 */
                                                  {0.400,300.*1E-4,0.9,0.7,0.7,1000.*1E-4,1000.*1E-4,-5000.*1E-8,0.95,0.,3.},  /*BIN 10*/ /* 10<pt<12 */
                                                  {0.400,300.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4, 10000.*1E-8,0.95,0.,3.},  /*BIN 11*/ /* 12<pt<16 */
                                                  {0.400,300.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4, 10000.*1E-8,0.90,0.,3.},  /*BIN 12*/ /* 16<pt<24 */
                                                  {0.400,300.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4, 10000.*1E-8,0.90,0.,3.}}; /*BIN 13*/ /* pt>24 */
  
//topomatic cut
  Float_t d0Topomatic[nptbins] = {2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.};
  RDHFD0Corr->Setd0MeasMinusExpCut(nptbins,d0Topomatic);
  
  //CREATE TRANSPOSE MATRIX...REVERSE INDICES as required by AliRDHFCuts
  Float_t **cutsMatrixTransposeStand=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++)cutsMatrixTransposeStand[iv]=new Float_t[nptbins];
 
  for (Int_t ibin=0;ibin<nptbins;ibin++){
    for (Int_t ivar = 0; ivar<nvars; ivar++){
      cutsMatrixTransposeStand[ivar][ibin]=cutsMatrixD0toKpiStand[ibin][ivar];      
    }
  }
 
  RDHFD0Corr->SetCuts(nvars,nptbins,cutsMatrixTransposeStand);
  RDHFD0Corr->SetUseSpecialCuts(kTRUE);
  RDHFD0Corr->SetRemoveDaughtersFromPrim(kTRUE);
  
  for(Int_t iv=0;iv<nvars;iv++) delete [] cutsMatrixTransposeStand[iv];
  delete [] cutsMatrixTransposeStand;
  cutsMatrixTransposeStand=NULL;
 
  //D0 pid settings
  Bool_t pidflag=kTRUE;
  RDHFD0Corr->SetUsePID(pidflag);
  if(pidflag) cout<<"PID is used"<<endl;
  else cout<<"PID is not used"<<endl;

  // PID SETTINGS
  AliAODPidHF* pidObj=new AliAODPidHF();
  //pidObj->SetName("pid4D0");
  Int_t mode=1;
  const Int_t nlims=2;
  Double_t plims[nlims]={0.6,0.8}; //TPC limits in momentum [GeV/c]
  Bool_t compat=kTRUE; //effective only for this mode
  Bool_t asym=kTRUE;
  Double_t sigmas[5]={2.,1.,0.,3.,0.}; //to be checked and to be modified with new implementation of setters by Rossella
  pidObj->SetAsym(asym);// if you want to use the asymmetric bands in TPC
  pidObj->SetMatch(mode);
  pidObj->SetPLimit(plims,nlims);
  pidObj->SetSigma(sigmas);
  pidObj->SetCompat(compat);
  pidObj->SetTPC(kTRUE);
  pidObj->SetTOF(kTRUE);
  pidObj->SetPCompatTOF(2.);
  pidObj->SetSigmaForTPCCompat(3.);
  pidObj->SetSigmaForTOFCompat(3.);
  pidObj->SetOldPid(kFALSE);
  RDHFD0Corr->SetPidHF(pidObj);
  RDHFD0Corr->SetUsePID(kTRUE);
  RDHFD0Corr->SetUseDefaultPID(kFALSE); //to use the AliAODPidHF
  RDHFD0Corr->SetLowPt(kFALSE);
  //RDHFD0Corr->SetMaximumPforPID(4.);

  //activate pileup rejection (for pp)
//  RDHFD0Corr->SetOptPileup(AliRDHFCuts::kRejectPileupEvent);

  TString cent="";
  //centrality selection
  if(isOnData) {
    RDHFD0Corr->SetMinCentrality(minc);
    RDHFD0Corr->SetMaxCentrality(maxc);
    cent=Form("%.0f%.0f",minc,maxc);
    RDHFD0Corr->SetUseCentrality(AliRDHFCuts::kCentZNA); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
  } else {
	RDHFD0Corr->SetUseCentrality(AliRDHFCuts::kCentOff);
	cent=Form("%.0f%.0f",minc,maxc);
  }

  //temporary
  //RDHFD0Corr->SetFixRefs();

  cout<<"This is the object I'm going to save:"<<endl;
  RDHFD0Corr->PrintAll();

  TFile* fout;
  if(isOnData) fout=new TFile("D0toKpiCuts_pPb2016_0_100_PreliminarySQM.root","recreate");   //set this!! 
  else fout=new TFile("D0toKpiCuts_pPb2016_0_100_PreliminarySQM_MC.root","recreate");   //set this!! 

  fout->cd();
  RDHFD0Corr->Write();
  fout->Close();

//____________________________________________________

  // Cuts for correlated tracks

  AliHFAssociatedTrackCuts* HFCorrelationCuts=new AliHFAssociatedTrackCuts();
  HFCorrelationCuts->SetName("AssociatedTrkCuts");
  HFCorrelationCuts->SetTitle("Cuts for associated track");
  Float_t eta = 0.8;

  // Set quality cuts on tracks
  AliESDtrackCuts *esdHadrCuts = new AliESDtrackCuts("AliESDHadrCuts","default");
  esdHadrCuts->SetRequireSigmaToVertex(kFALSE);
  esdHadrCuts->SetRequireTPCRefit(kTRUE);
  esdHadrCuts->SetRequireITSRefit(kFALSE);
  esdHadrCuts->SetMinNClustersITS(2);
//esdHadrCuts->SetMinNClustersTPC(70);
  esdHadrCuts->SetMinNCrossedRowsTPC(70);
  esdHadrCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  esdHadrCuts->SetMaxChi2PerClusterTPC(4);
//esdHadrCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  esdHadrCuts->SetMaxDCAToVertexZ(1.);
  esdHadrCuts->SetMaxDCAToVertexXY(1.);
  esdHadrCuts->SetEtaRange(-eta,eta);
  esdHadrCuts->SetPtRange(0.3,1.e10);

  HFCorrelationCuts->AddTrackCuts(esdHadrCuts);

  // Set kinematics cuts for AOD track 
  const int nofcuts = 4;
  Float_t* trackcutsarray;
  trackcutsarray=new Float_t[nofcuts];
  trackcutsarray[0] = 0.3;//track min pt
  trackcutsarray[1] = 10000.;//track max pt
  trackcutsarray[2] = -1000000000.;//track min impact parameter. DON'T put 0 since default value is -999999. and it would skip all tracks if d0 is not calculated!
  trackcutsarray[3] = 10000.;//track max impact parameter
  HFCorrelationCuts->SetNVarsTrack(nofcuts);
  HFCorrelationCuts->SetAODTrackCuts(trackcutsarray);

  HFCorrelationCuts->SetCharge(0); // -1/+1 to look for opposite/same charge, 0 no charge selection 
  HFCorrelationCuts->SetFilterBit(AliAODTrack::kTrkTPCOnly); // set 0 for analysis with AOD from 2010, il 4 è kTrkGlobalNoDCA
	
  // Set kinematics cuts for AOD v0 
  const int nofcuts2 = 7;
  Float_t* vzerocutsarray;
  vzerocutsarray=new Float_t[nofcuts2];
  vzerocutsarray[0] = 0.2; // max dca between two daugters (cm)
  vzerocutsarray[1] = 2; //  max chi square
  vzerocutsarray[2] = 2.; // min decay length (cm) 
  vzerocutsarray[3] = 1000; // max decay length (cm)
  vzerocutsarray[4] = 1000.; // max opening angle between two daugters
  vzerocutsarray[5] = 0; // min pt of k0 (GeV/c)
  vzerocutsarray[6] = 0.9; // set eta acceptance
  HFCorrelationCuts->SetNVarsVzero(nofcuts2);
  HFCorrelationCuts->SetAODvZeroCuts(vzerocutsarray);
	
  // Set PID
  Int_t mode =1;
  Double_t ptlimit[2] = {0.6,0.8};
  AliAODPidHF* pidObj=new AliAODPidHF();
  pidObj->SetMatch(1);  //A.Rossi mode
  pidObj->SetAsym(kTRUE);
  pidObj->SetPLimit(ptlimit,2);
  pidObj->SetSigma(0,2.);  //TPC sigma, in three pT ranges
  pidObj->SetSigma(1,1.);
  pidObj->SetSigma(2,0.);  
  pidObj->SetSigma(3,3.);  //TOF sigma, whole pT range
  pidObj->SetPCompatTOF(1.5);
  pidObj->SetSigmaForTPCCompat(3.);
  pidObj->SetSigmaForTOFCompat(3.);
  pidObj->SetTPC(kTRUE);
  pidObj->SetTOF(kTRUE);
  pidObj->SetCompat(kTRUE);
  HFCorrelationCuts->SetPidHF(pidObj);

  //Event Mixing settings
  HFCorrelationCuts->SetMaxNEventsInPool(200);
  HFCorrelationCuts->SetMinNTracksInPool(10000); //increase to 5000? but modifying in HFCorrelator fTargetDepth
  HFCorrelationCuts->SetMinEventsToMix(1);  //reduce to 1?
  HFCorrelationCuts->SetTargetFracTracks(0.0025);
  HFCorrelationCuts->SetNofPoolBins(3,3);

  Double_t MBins[]={0,35,55,500};
  Double_t * MultiplicityBins = MBins;
  Double_t ZBins[]={-10,-1.5,3.5,10};
  Double_t *ZVrtxBins = ZBins;

  HFCorrelationCuts->SetPoolBins(ZVrtxBins,MultiplicityBins);

  //Select MC production process
  HFCorrelationCuts->SetNofMCEventTypes(3); //1 for only one
  Int_t MCEvType[]={28,53,68}; //28 = Flavour excitation, 53 = Pair production, 68 = Gluon splitting
  Int_t *MCEvTypeArray = MCEvType;
  HFCorrelationCuts->SetMCEventTypes(MCEvTypeArray);

  //Define sideband edges (if set externally in the AddTask file)
  Double_t LSBLow[14] = {0.,0.,1.7928,1.7928,1.7768,1.7728,1.7648,1.7488,1.7448,1.7428,1.7428,1.7048,1.7048,1.7048}; //to be filled looking at results from invariant mass fits!
  Double_t LSBUpp[14] = {0.,0.,1.8288,1.8288,1.8208,1.8208,1.8128,1.8088,1.8048,1.8048,1.8048,1.7568,1.7568,1.7568};
  Double_t RSBLow[14] = {0.,0.,1.9008,1.9008,1.9088,1.9128,1.9168,1.9288,1.9288,1.9288,1.9288,1.9728,1.9728,1.9728};
  Double_t RSBUpp[14] = {0.,0.,1.9408,1.9408,1.9528,1.9608,1.9688,1.9848,1.9888,1.9928,1.9928,2.0768,2.0768,2.0768}; 

  TVectorD vLSBLow(14,LSBLow);
  TVectorD vLSBUpp(14,LSBUpp);
  TVectorD vRSBLow(14,RSBLow);
  TVectorD vRSBUpp(14,RSBUpp);

  // Save to *.root file
  HFCorrelationCuts->PrintAll();
  TFile* fout=new TFile("AssocPartCuts_Std_pPb2016_0_100_PreliminarySQM.root","recreate");   //set this!! 
  fout->cd();
  HFCorrelationCuts->Write();
  vLSBLow.Write("vLSBLow");
  vLSBUpp.Write("vLSBUpp");
  vRSBLow.Write("vRSBLow");
  vRSBUpp.Write("vRSBUpp");
  fout->Close();

}

