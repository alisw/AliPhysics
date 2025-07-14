//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Konstantin.Mikhyalov@cern.ch
// Jan 2025 --> All wagon before 2025 failed
//              ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// pPb5.02(run2) -> K+K+ & K-K-
// PID is the same as in PbPb2.76->K+K- 014 wagon !!!
// 0.2-0.5,0.5-1.5; 0-20,20-40,40-90%
// event shape: St:
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#if !defined(__CINT__) || defined(__MAKECINT_)

#include "AliFemtoManager.h"
//#include "AliFemtoEventReaderAOD.h"
#include "AliFemtoEventReaderAODMultSelection.h"
#include "AliFemtoVertexMultAnalysis.h"
//#include "AliFemtoBasicEventCut.h" //w/o sphericity
#include "AliFemtoSphericityEventCut.h" //w/ sphericity
#include "AliFemtoKpm45TrackCut.h"
#include "AliFemtoCutMonitorParticlePID.h"
#include "AliFemtoPairCutRadialDistanceKKdist.h"
#include "AliFemtoQinvCorrFctn.h"
#include "AliFemtoKTPairCut.h"
#include "AliFemtoTPCInnerCorrFctn.h"

#endif

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis() {

  const double m_pi = 0.13956995;
  const double m_K  = 0.493677;
  const int cMu=3; //numbwr of multiplicity bins
  int runmults[cMu] = {1, 1, 1};
  double multbins[cMu+1] = {0.01, 200, 400, 900};//LEGO Train
  //double multbins[cMu+1] = {0, 50, 150, 900};// 12 AOD.root file
  //double multbins[cMu+1] = {0, 50, 100, 200};// 1 AOD.root file
  const char *chrgs[2] = { "Kp", "Km"};
  const char *spart[2] = {"Kp", "Km"};
  const char *sKK[2] = {"KpKp", "KmKm"};//CF names
  const int zar[2]={+1,-1};

  //const double S_T_MIMA[2]={0.0,1.0};//all events
  //const double S_T_MIMA[2]={0.0,0.3};//jetty events
  const double S_T_MIMA[2]={0.7,1.0};//sperical events
  const double DCA[2] = {0.135,0.130};//xy,z in mm

  const Bool_t cf_kT = kTRUE;//CF for a few k_T bins
  const Int_t  cKt=2;
  //double ktrng[cKt+1] = {0.2, 0.5, 1.0, 1.5};//orig 0.2-0.5,0.5-1.0
  double ktrng[cKt+1] = {0.15, 0.5, 1.2};//like Lyuda's

  int run3d = -1/*1*/;
  int runshlcms = 0;
  double shqmax;
  int nbinssh = 200;
  //int nbinssh = 100;

  //if (runshlcms) shqmax = 2.0;
  if (runshlcms) shqmax = 0.25;
  else shqmax = 2.0;

  //pPb 5.02 config reader:
  //Both AliFemtoEventReaderAODMultSelection *Reader and AliFemtoEventReaderAOD *Reader are the same
  AliFemtoEventReaderAOD *Reader = new AliFemtoEventReaderAODMultSelection();
  Reader->SetFilterBit(8); //8 - hybrid (TPC+ITS?), 7 - TPC. PbPb: AOD FilterBit 128 which means that TPC only tracks constrained to SPD primary vertex are used
    Reader->SetEPVZERO(kTRUE); // How do we use EventPlaneEngle???
    Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);//no difference if Reader->SetUseMultiplicity(AliFemtoEventReaderAODMultSelection::kCentralityV0A);//Lena
    //Reader->SetUseMultiplicity(AliFemtoEventReaderAODMultSelection::kCentralityV0A);//Lena
    Reader->SetCentralityFlattening(kFALSE/*kTRUE*/);// This should be checked!!
    Reader->SetDCAglobalTrack(kTRUE);
    Reader->SetReadV0(0);//0 is OK, but 1 *** Break *** segmentation violation. It was 1 in Lena's config. Did Lena study V0? Lena took this from some of Polish configuration without any reasons.
   
  AliFemtoManager* Manager=new AliFemtoManager();
    Manager->SetEventReader(Reader);
    
  const int iz=20;//size of arrays
  AliFemtoVertexMultAnalysis    *mixi[iz];
  AliFemtoSphericityEventCut    *bcut[iz]; //w/ spher:
  //AliFemtoBasicEventCut         *bcut[iz];//w/o sphericity
  
  AliFemtoCutMonitorEventSphericity *smon[2][iz];//sphericity monitor
  AliFemtoCutMonitorEventVertex     *vmon[2][iz];//vertex monitor
  AliFemtoCutMonitorEventMult       *mul_[2][iz];//multiplicity
  AliFemtoKpm45TrackCut             *pid_[iz];// pid_[2][iz] in K+K- case
  AliFemtoPairCutRadialDistanceKKdist  *ttc[iz];
  AliFemtoQinvCorrFctn                 *cqi[iz*cKt];//1d CF(Qinv): all kT + cKt bins

  AliFemtoCutMonitorParticleYPt *mdca[2][iz];//[pass-fail][iz]
  AliFemtoCutMonitorParticlePID *mpid[2][iz];//[pass-fail][iz]
  
  AliFemtoKTPairCut             *pcut[iz*cKt];//Pair cut like k_T range
  AliFemtoTPCInnerCorrFctn      *mdis[iz*cKt];//distance at the entrance to the TPC between two tracks
  
  int ix = 0;
  for (int m=0; m<cMu; m++) {
    if (runmults[m] == 0) continue;
    for (int zz=0;zz<2;zz++) {//indeks znaka zaryada dlya tozhdestvennyh
    //just print
    for(int d=0;d<5;d++) cout<<"====================\n";
    cout<<"====================Multiplisity="<<m<<"\n";
    ix = m+zz*3;//ichg*3+imult;???

      mixi[ix] = new AliFemtoVertexMultAnalysis(10, -10.0, 10.0, 4, multbins[m], multbins[m+1]);
      mixi[ix]->SetNumEventsToMix(5);
      mixi[ix]->SetMinSizePartCollection(1);
      mixi[ix]->SetVerboseMode(kTRUE/*kFALSE*/);
	
      bcut[ix] = new AliFemtoSphericityEventCut();//bcut[ix] = new AliFemtoBasicEventCut();//spher: 
      bcut[ix]->SetEventMult(0,10000);
      bcut[ix]->SetVertZPos(-10,10);
      bcut[ix]->SetStMin(S_T_MIMA[0]);
      bcut[ix]->SetStMax(S_T_MIMA[1]);
      
      smon[0][ix] = new AliFemtoCutMonitorEventSphericity(Form("_Pass_M%i",m));
      smon[1][ix] = new AliFemtoCutMonitorEventSphericity(Form("_Fail_M%i",m));
      bcut[ix]->AddCutMonitor(smon[0][ix],smon[1][ix]);

	
      //AliFemtoCutMonitorEventVertex class gives empty XY and wide Z vertex hists ???
      //vmon[0][ix] = new AliFemtoCutMonitorEventVertex(Form("_Pass_M%i",m));
      //vmon[1][ix] = new AliFemtoCutMonitorEventVertex(Form("_Fail_M%i",m));
      //bcut[ix]->AddCutMonitor(vmon[0][ix],vmon[1][ix]);

      mul_[0][ix] = new AliFemtoCutMonitorEventMult(Form("_Pass_M%i",m)); 
      mul_[1][ix] = new AliFemtoCutMonitorEventMult(Form("_Fail_M%i",m));
      bcut[ix]->AddCutMonitor(mul_[0][ix],mul_[1][ix]);
      //identical particles::no loop over charge
      //for(int pa=0;pa<1;pa++) {
	pid_[ix] = new AliFemtoKpm45TrackCut();
	pid_[ix]->SetCharge(zar[zz]);//zar[0]=+1 zar[1]=-1 
	//The low limit of TOF using	
	//pid_[ix]->SetPMinTOF45(0.405); //Since March6,2025

	pid_[ix]->SetPt(0.14,1.5);
	pid_[ix]->SetEta(-0.8,0.8);
	pid_[ix]->SetMass(m_K);
	pid_[ix]->SetMostProbableKaon();
	pid_[ix]->SetNsigmaTPCle250(2.0);    //K+K- standard PID -->
	pid_[ix]->SetNsigmaTPC250_400(2.0);  //
	pid_[ix]->SetNsigmaTPC400_450(1.0);  //to remove e-peak decrease from 3.0(default) to smaller value = 1.0(orig!!!!)
	pid_[ix]->SetNsigmaTPC450_500(2.0);  //orig->2.0 --> does not use !!!
	pid_[ix]->SetNsigmaTPCge500(3.0);    //orig->3.0
	pid_[ix]->SetNsigmaTOF500_800(2.0);  //
	pid_[ix]->SetNsigmaTOF800_1000(1.5); //
	pid_[ix]->SetNsigmaTOFge1000(1.0);   //K+K- standard PID --<
	
	pid_[ix]->SetminTPCncls(80);         //K+K- track quality cuts -->
	pid_[ix]->SetRemoveKinks(kTRUE);     //
	pid_[ix]->SetMaxTPCChiNdof(4.0);     //
      //pid_[ix]->SetMaxITSChiNdof(36);//	 //  
	pid_[ix]->SetLabel(kFALSE);          //
	pid_[ix]->SetMaxImpactXY(DCA[0]);    //
	pid_[ix]->SetMaxImpactZ( DCA[1]);    //K+K- track quality cuts --<

	mdca[0][ix] = new AliFemtoCutMonitorParticleYPt(Form("_Pass_%s_M%i",spart[zz],m),m_K);
	mdca[1][ix] = new AliFemtoCutMonitorParticleYPt(Form("_Fail_%s_M%i",spart[zz],m),m_K);
	pid_[ix]->AddCutMonitor(mdca[0][ix],mdca[1][ix]);

	mpid[0][ix] = new AliFemtoCutMonitorParticlePID(Form("_Pass_%s_M%i", spart[zz],m),1);//1 - kaon time
	mpid[1][ix]  = new AliFemtoCutMonitorParticlePID(Form("_Fail_%s_M%i",spart[zz],m),1);
	pid_[ix]->AddCutMonitor(mpid[0][ix],mpid[1][ix]);

	//}//for(int pa=0;pa<1;pa++) {
      ttc[ix] = new AliFemtoPairCutRadialDistanceKKdist(); //K+K- two track cuts -->
      ttc[ix]->SetShareQualityMax(1.0);                    //
      ttc[ix]->SetShareFractionMax(0.05);                  //
      ttc[ix]->SetAverageSeparation(3.0);                  //
      ttc[ix]->SetRemoveSameLabel(kFALSE);                 //K+K- two track cuts --<
	  
	  
      mixi[ix]->SetEventCut(bcut[ix]);
       mixi[ix]->SetFirstParticleCut(pid_[ix]);//The same since pid the same
      mixi[ix]->SetSecondParticleCut(pid_[ix]);//for identical particles
      mixi[ix]->SetPairCut(ttc[ix]);
	  
      //CF vs Qinv (without kT bins)
      cqi[ix] = new AliFemtoQinvCorrFctn(Form("C_%s_Mu_%i", sKK[zz], m),nbinssh,0.0,shqmax);
      mixi[ix]->AddCorrFctn(cqi[ix]);
      
      //CF vs Qinv for k_T
      int ik=0;
      for(int k=0;k<cKt;k++){
	if(!cf_kT) continue;
	ik=ix*cKt+k+cMu;//"+1" to be different with cqi[ix]
	cout<<"------------------> ik= "<<ik<<"  ix=m+3*zz="<<ix<<endl;
	pcut[ik]=new AliFemtoKTPairCut(ktrng[k],ktrng[k+1]);
	cqi[ik] =new AliFemtoQinvCorrFctn(Form("C_%s_Mu_%i_kT_%i", sKK[zz], m,k),nbinssh,0.0,shqmax);
	cqi[ik]->SetPairSelectionCut(pcut[ik]);
	mixi[ix]->AddCorrFctn(cqi[ik]);
	
	mdis[ik]=new AliFemtoTPCInnerCorrFctn(Form("_DarR_Mu_%i_kT_%i",m,k),nbinssh,0.0,shqmax);
	mdis[ik]->SetRadius(1.2);//DeltaEta-DeltaPhi* calculated at this Radius
	mixi[ix]->AddCorrFctn(mdis[ik]);
	/*
	  shar[ktm] = new AliFemtoShareQualityCorrFctn(Form("_Share_%s_M_%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);//A correlation function that saves the     ///
/// amount of sharing and splitting hits per pair as a function of qinv 
	  shar[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
	  mixi[ik]->AddCorrFctn(shar[ik]]);
	 */
	
      }
      Manager->AddAnalysis(mixi[ix]);
      
    }//for (int zz=0;zz<2;zz++) {//indeks znaka zaryada dlya tozhdestvennyh
  }//  for (int m=0; m<cMu; m++) {

  return Manager;
}                         
                      
