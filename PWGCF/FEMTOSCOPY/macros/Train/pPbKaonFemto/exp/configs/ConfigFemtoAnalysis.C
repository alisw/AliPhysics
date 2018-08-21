//Tomasz Lehmann

#if !defined(__CINT__) || defined(__MAKECINT_)
#include "AliFemtoManager.h"
#include "AliFemtoEventReaderESDChain.h"
#include "AliFemtoEventReaderESDChainKine.h"
#include "AliFemtoEventReaderAODChain.h"
#include "AliFemtoSimpleAnalysis.h"
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoESDTrackCut.h"
#include "AliFemtoKKTrackCut.h"
#include "AliFemtoCorrFctn.h"
#include "AliFemtoCutMonitorParticleYPt.h"
#include "AliFemtoCutMonitorParticleVertPos.h"
#include "AliFemtoCutMonitorParticleMomRes.h"
#include "AliFemtoCutMonitorParticlePID.h"
#include "AliFemtoCutMonitorEventMult.h"
#include "AliFemtoCutMonitorEventVertex.h"
#include "AliFemtoShareQualityTPCEntranceSepPairCut.h"
#include "AliFemtoPairCutAntiGamma.h"
#include "AliFemtoPairCutRadialDistance.h"
#include "AliFemtoQinvCorrFctn.h"
#include "AliFemtoShareQualityCorrFctn.h"
#include "AliFemtoTPCInnerCorrFctn.h"
#include "AliFemtoVertexMultAnalysis.h"
#include "AliFemtoCorrFctn3DSpherical.h"
#include "AliFemtoChi2CorrFctn.h"
#include "AliFemtoCorrFctnTPCNcls.h"
#include "AliFemtoBPLCMS3DCorrFctn.h"
#include "AliFemtoCorrFctn3DLCMSSym.h"
#include "AliFemtoModelBPLCMSCorrFctn.h"
#include "AliFemtoModelCorrFctn3DSpherical.h"
#include "AliFemtoModelGausLCMSFreezeOutGenerator.h"
#include "AliFemtoModelGausRinvFreezeOutGenerator.h"
#include "AliFemtoModelManager.h"
#include "AliFemtoModelWeightGeneratorBasic.h"
#include "AliFemtoModelWeightGeneratorLednicky.h"
#include "AliFemtoCorrFctnDirectYlm.h"
#include "AliFemtoModelCorrFctnDirectYlm.h"
#include "AliFemtoModelCorrFctnSource.h"
#include "AliFemtoCutMonitorParticlePtPDG.h"
#include "AliFemtoKTPairCut.h"
#include "AliFemtoCutMonitorCollections.h"
#endif

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis() {

  double PionMass = 0.13956995;
  double KaonMass = 0.493677;
	
   int runmults[3] = {1, 1, 1};
    int multbins[4] = {0.01, 200, 400, 900};

  int runch[2] = {1, 1}; // Why?
  const char *chrgs[2] = { "Kp", "Km"};

  int runktdep = 1;
  double ktrng[3] = {0.2, 0.5, 1.0};

  int run3d = 0; // Do 3D cartesian analysis?
  //int runshlcms = 1;
  int runshlcms = 0;
  double shqmax;
  //int nbinssh = 200;
  int nbinssh = 100;

  //if (runshlcms) shqmax = 2.0;
  if (runshlcms) shqmax = 0.25;
  else shqmax = 2.0;

  AliFemtoEventReaderAODChain *Reader = new AliFemtoEventReaderAODChain();
  Reader->SetUseMultiplicity(AliFemtoEventReaderAODChain::kCentrality);
  Reader->SetFilterBit(5);
  //Reader->SetDCAglobalTrack(kTRUE);
  Reader->SetpA2013(kTRUE);
  
  AliFemtoManager* Manager=new AliFemtoManager();
  Manager->SetEventReader(Reader);

  AliFemtoVertexMultAnalysis    *anetaphitpc[20];
  AliFemtoVertexMultAnalysis    *anetaphitpc1[20];
  AliFemtoVertexMultAnalysis    *anetaphitpc2[20];
  AliFemtoBasicEventCut         *mecetaphitpc[20];
  AliFemtoCutMonitorEventMult   *cutPassEvMetaphitpc[20];
  AliFemtoCutMonitorEventMult   *cutFailEvMetaphitpc[20];
  AliFemtoCutMonitorEventVertex *cutPassEvVetaphitpc[20];
  AliFemtoCutMonitorEventVertex *cutFailEvVetaphitpc[20];
  AliFemtoCutMonitorCollections   *cutPassColletaphitpc[20];
  AliFemtoCutMonitorCollections   *cutFailColletaphitpc[20];
  AliFemtoKKTrackCut           *dtc1etaphitpc[20];
  AliFemtoKKTrackCut           *dtc2etaphitpc[20];
 // AliFemtoKKTrackCut           *dtc1etaphitpcstrongcuts[20]; //dodane
 // AliFemtoKKTrackCut           *dtc2etaphitpcstrongcuts[20]; //dodane
  AliFemtoCutMonitorParticleYPt *cutPass1YPtetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutFail1YPtetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutPass1PIDetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutFail1PIDetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutPass2YPtetaphitpc[20];
  AliFemtoCutMonitorParticleYPt *cutFail2YPtetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutPass2PIDetaphitpc[20];
  AliFemtoCutMonitorParticlePID *cutFail2PIDetaphitpc[20];
  AliFemtoPairCutAntiGamma      *sqpcetaphitpc[20];
//  AliFemtoPairCutAntiGamma      *sqpcetaphitpc1[20]; //dodane
  AliFemtoCorrFctnDirectYlm     *cylmetaphitpc[20];
  AliFemtoCorrFctnDEtaDPhi      *cdedpetaphi[20];
  AliFemtoChi2CorrFctn          *cchiqinvetaphitpc[20];
  AliFemtoCorrFctnGammaMonitor  *cgamma[20*10];
  AliFemtoKTPairCut             *ktpcuts[20*8];
  AliFemtoCorrFctnDirectYlm     *cylmkttpc[20*8];
  AliFemtoQinvCorrFctn          *cqinvkttpc[20*8];
  AliFemtoCorrFctn3DLCMSSym     *cq3dlcmskttpc[20*8];
  AliFemtoCorrFctnTPCNcls       *cqinvnclstpc[20];
  AliFemtoShareQualityCorrFctn  *cqinvsqtpc[20*10];
  AliFemtoChi2CorrFctn          *cqinvchi2tpc[20];
  AliFemtoTPCInnerCorrFctn      *cqinvinnertpc[20*10];

  // *** Begin analysis ***
  int aniter = 0;
  int aniter1=0;
  int aniter2=0;
  int ichg=0;
  for (int imult=0; imult<3; imult++) {
    if (runmults[imult]) {
        if (runch[ichg]) {
          aniter = ichg*3+imult;

          anetaphitpc[aniter] = new AliFemtoVertexMultAnalysis(10, -10.0, 10.0, 4, multbins[imult], multbins[imult+1]);
          anetaphitpc[aniter]->SetNumEventsToMix(5);
          anetaphitpc[aniter]->SetMinSizePartCollection(1);
          anetaphitpc[aniter]->SetVerboseMode(kTRUE); 

          mecetaphitpc[aniter] = new AliFemtoBasicEventCut();
          mecetaphitpc[aniter]->SetEventMult(0,10000);
          mecetaphitpc[aniter]->SetVertZPos(-10,10);

          cutPassEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cut1Pass%stpcM%i", chrgs[ichg], imult));
          cutFailEvMetaphitpc[aniter] = new AliFemtoCutMonitorEventMult(Form("cut1Fail%stpcM%i", chrgs[ichg], imult));
          mecetaphitpc[aniter]->AddCutMonitor(cutPassEvMetaphitpc[aniter], cutFailEvMetaphitpc[aniter]);

          cutPassEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cut1Pass%stpcM%i", chrgs[ichg], imult));
          cutFailEvVetaphitpc[aniter] = new AliFemtoCutMonitorEventVertex(Form("cut1Fail%stpcM%i", chrgs[ichg], imult));
          mecetaphitpc[aniter]->AddCutMonitor(cutPassEvVetaphitpc[aniter], cutFailEvVetaphitpc[aniter]);

          cutPassColletaphitpc[aniter] = new AliFemtoCutMonitorCollections(Form("cut1Pass%stpcM%i", chrgs[ichg], imult));
              cutFailColletaphitpc[aniter] = new AliFemtoCutMonitorCollections(Form("cut1Fail%stpcM%i", chrgs[ichg], imult));
              mecetaphitpc[aniter]->AddCutMonitor(cutPassColletaphitpc[aniter], cutFailColletaphitpc[aniter]);
            
            
//-----------------------1 particle-------------------------------------------<
            dtc1etaphitpc[aniter] = new AliFemtoKKTrackCut();
            dtc1etaphitpc[aniter]->SetCharge(1.0);
            dtc1etaphitpc[aniter]->SetPt(0.14,1.5);
            dtc1etaphitpc[aniter]->SetEta(-0.8,0.8);
          //PID method
            dtc1etaphitpc[aniter]->SetMass(KaonMass);
            dtc1etaphitpc[aniter]->SetMostProbableKaon();

         dtc1etaphitpc[aniter]->SetNsigmaTPCle250(2.0);
         dtc1etaphitpc[aniter]->SetNsigmaTPC250_400(2.0);
          dtc1etaphitpc[aniter]->SetNsigmaTPC400_450(2.0);
          dtc1etaphitpc[aniter]->SetNsigmaTPC450_500(2.0);
          dtc1etaphitpc[aniter]->SetNsigmaTPCge500(3.0); 
          dtc1etaphitpc[aniter]->SetNsigmaTOF500_800(2.0);
          dtc1etaphitpc[aniter]->SetNsigmaTOF800_1000(1.5);
          dtc1etaphitpc[aniter]->SetNsigmaTOFge1000(1.0);
          //------------------- November 2013 ----------------------------------->
        
          dtc1etaphitpc[aniter]->SetRemoveKinks(kTRUE);
  
          dtc1etaphitpc[aniter]->SetLabel(kFALSE);
    //----------------------2particle----------------------< KR 
            // dtc2etaphitpc[aniter] = new AliFemtoESDTrackCut();
            dtc2etaphitpc[aniter]=new AliFemtoKKTrackCut();
            dtc2etaphitpc[aniter]->SetCharge(-1.0);
            dtc2etaphitpc[aniter]->SetPt(0.14,1.5);
            dtc2etaphitpc[aniter]->SetEta(-0.8,0.8);
          //PID method
            dtc2etaphitpc[aniter]->SetMass(KaonMass);
            dtc2etaphitpc[aniter]->SetMostProbableKaon();
          //dtc2etaphitpc[aniter]->SetPIDMethod(AliFemtoESDTrackCut::kContour);
            //------------------- November 2013 -----------------------------------< 
          // new cuts to remove electron (do not take into analysis if 400<p<500) 
         dtc2etaphitpc[aniter]->SetNsigmaTPCle250(2.0);
         dtc2etaphitpc[aniter]->SetNsigmaTPC250_400(2.0);
          dtc2etaphitpc[aniter]->SetNsigmaTPC400_450(2.0);
          dtc2etaphitpc[aniter]->SetNsigmaTPC450_500(2.0);
          dtc2etaphitpc[aniter]->SetNsigmaTPCge500(3.0);    
          // new cuts are stronger, better separation of pion in TOF 
          // when momentum is greater then 800 MeV/c
          dtc2etaphitpc[aniter]->SetNsigmaTOF500_800(2.0);
          dtc2etaphitpc[aniter]->SetNsigmaTOF800_1000(1.5);
          dtc2etaphitpc[aniter]->SetNsigmaTOFge1000(1.0);
          //------------------- November 2013 ----------------------------------->
          //Track quality cuts
      
          dtc2etaphitpc[aniter]->SetRemoveKinks(kTRUE);
          	  
          dtc2etaphitpc[aniter]->SetLabel(kFALSE);
          cutPass1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cut1Pass1%stpcM%i", chrgs[ichg], imult), 0.493677);
          cutFail1YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cut1Fail1%stpcM%i", chrgs[ichg], imult), 0.493677);
            dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1YPtetaphitpc[aniter], cutFail1YPtetaphitpc[aniter]);
            
            cutPass2YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cut1Pass2%stpcM%i", chrgs[ichg+1], imult), 0.493677); //ichg+1 --> ichg=1 --> charge = -1
            cutFail2YPtetaphitpc[aniter] = new AliFemtoCutMonitorParticleYPt(Form("cut1Fail2%stpcM%i", chrgs[ichg+1], imult), 0.493677);
          
          dtc2etaphitpc[aniter]->AddCutMonitor(cutPass2YPtetaphitpc[aniter], cutFail2YPtetaphitpc[aniter]);
/*****************************************************/

          cutPass1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cut1Pass1%stpcM%i", chrgs[ichg], imult),1);
          cutFail1PIDetaphitpc[aniter] = new AliFemtoCutMonitorParticlePID(Form("cut1Fail1%stpcM%i", chrgs[ichg], imult),1);
          dtc1etaphitpc[aniter]->AddCutMonitor(cutPass1PIDetaphitpc[aniter], cutFail1PIDetaphitpc[aniter]);
 
 //cuts on gamma------------------------------------------------         
          sqpcetaphitpc[aniter] = new AliFemtoPairCutAntiGamma();
          sqpcetaphitpc[aniter]->SetRemoveSameLabel(kFALSE);
          sqpcetaphitpc[aniter]->SetMaxEEMinv(0.001); 
          sqpcetaphitpc[aniter]->SetMaxThetaDiff(0.02);
          sqpcetaphitpc[aniter]->SetTPCEntranceSepMinimum(0); 
 //-------------------------------------------------------------
       
/*****************************************************/
          anetaphitpc[aniter]->SetEventCut(mecetaphitpc[aniter]);
          anetaphitpc[aniter]->SetFirstParticleCut(dtc1etaphitpc[aniter]);
          anetaphitpc[aniter]->SetSecondParticleCut(dtc2etaphitpc[aniter]); //druga czastka
          anetaphitpc[aniter]->SetPairCut(sqpcetaphitpc[aniter]); 
/*****************************************************/
          //Qinv (without kT bins)
          cqinvkttpc[aniter] = new AliFemtoQinvCorrFctn(Form("1cqinv%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
          anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[aniter]);

          //3D cartesian (without kT bins)
          if(run3d){
            cq3dlcmskttpc[aniter] = new AliFemtoCorrFctn3DLCMSSym(Form("1cq3d%stpcM%i", chrgs[ichg], imult),100,0.5);
            anetaphitpc[aniter]->AddCorrFctn(cq3dlcmskttpc[aniter]);
          }

/*****************************************************/
          if (runktdep) {
            int ktm;
            for (int ikt=0; ikt<2; ikt++) {
              ktm = aniter*2 + ikt;
              ktpcuts[ktm] = new AliFemtoKTPairCut(ktrng[ikt], ktrng[ikt+1]);
/*****************************************************/
              cqinvkttpc[ktm] = new AliFemtoQinvCorrFctn(Form("1cqinv%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0, shqmax);
              cqinvkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
              anetaphitpc[aniter]->AddCorrFctn(cqinvkttpc[ktm]);

              cqinvsqtpc[ktm] = new AliFemtoShareQualityCorrFctn(Form("1cqinvsq%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
              cqinvsqtpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
              anetaphitpc[aniter]->AddCorrFctn(cqinvsqtpc[ktm]);
/*****************************************************/
              cqinvinnertpc[ktm] = new AliFemtoTPCInnerCorrFctn(Form("1cqinvinner%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
              cqinvinnertpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
              cqinvinnertpc[ktm]->SetRadius(1.6);
              anetaphitpc[aniter]->AddCorrFctn(cqinvinnertpc[ktm]);
              cgamma[aniter] = new AliFemtoCorrFctnGammaMonitor(Form("1cgammaM%ikT%i", imult, ikt),200,200);
              anetaphitpc[aniter]->AddCorrFctn(cgamma[aniter]);

              if (run3d) {
            cq3dlcmskttpc[ktm] = new AliFemtoCorrFctn3DLCMSSym(Form("1cq3d%stpcM%ikT%i", chrgs[ichg], imult, ikt),60,0.5);
            cq3dlcmskttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
            anetaphitpc[aniter]->AddCorrFctn(cq3dlcmskttpc[ktm]);
              }
            }
          }

          cdedpetaphi[aniter] = new AliFemtoCorrFctnDEtaDPhi(Form("1cdedp%stpcM%i", chrgs[ichg], imult),39, 39);
          anetaphitpc[aniter]->AddCorrFctn(cdedpetaphi[aniter]);

          Manager->AddAnalysis(anetaphitpc[aniter]);	
        }
    }
  }
 //===========================================strong cuts analysis===============================================================
    for (int imult=0; imult<3; imult++) {
    if (runmults[imult]) {
        if (runch[ichg]) {
          aniter1 = ichg*3+imult;

          anetaphitpc1[aniter1] = new AliFemtoVertexMultAnalysis(10, -10.0, 10.0, 4, multbins[imult], multbins[imult+1]);
          anetaphitpc1[aniter1]->SetNumEventsToMix(5);
          anetaphitpc1[aniter1]->SetMinSizePartCollection(1);
          anetaphitpc1[aniter1]->SetVerboseMode(kTRUE); 

          mecetaphitpc[aniter1] = new AliFemtoBasicEventCut();
          mecetaphitpc[aniter1]->SetEventMult(0,10000);
          mecetaphitpc[aniter1]->SetVertZPos(-10,10);

          cutPassEvMetaphitpc[aniter1] = new AliFemtoCutMonitorEventMult(Form("cut2Pass%stpcM%i", chrgs[ichg], imult));
          cutFailEvMetaphitpc[aniter1] = new AliFemtoCutMonitorEventMult(Form("cut2Fail%stpcM%i", chrgs[ichg], imult));
          mecetaphitpc[aniter1]->AddCutMonitor(cutPassEvMetaphitpc[aniter1], cutFailEvMetaphitpc[aniter1]);

          cutPassEvVetaphitpc[aniter1] = new AliFemtoCutMonitorEventVertex(Form("cut2Pass%stpcM%i", chrgs[ichg], imult));
          cutFailEvVetaphitpc[aniter1] = new AliFemtoCutMonitorEventVertex(Form("cut2Fail%stpcM%i", chrgs[ichg], imult));
          mecetaphitpc[aniter1]->AddCutMonitor(cutPassEvVetaphitpc[aniter1], cutFailEvVetaphitpc[aniter1]);

          cutPassColletaphitpc[aniter1] = new AliFemtoCutMonitorCollections(Form("cut2Pass%stpcM%i", chrgs[ichg], imult));
              cutFailColletaphitpc[aniter1] = new AliFemtoCutMonitorCollections(Form("cut2Fail%stpcM%i", chrgs[ichg], imult));
              mecetaphitpc[aniter1]->AddCutMonitor(cutPassColletaphitpc[aniter1], cutFailColletaphitpc[aniter1]);
            
//strong cuts-------------------------------------------------------------------------
  dtc1etaphitpc[aniter1] = new AliFemtoKKTrackCut();
            dtc1etaphitpc[aniter1]->SetCharge(1.0);
            dtc1etaphitpc[aniter1]->SetPt(0.14,0.45);
            dtc1etaphitpc[aniter1]->SetEta(-0.8,0.8);
            dtc1etaphitpc[aniter1]->SetMass(KaonMass);
            dtc1etaphitpc[aniter1]->SetMostProbableKaon();

         dtc1etaphitpc[aniter1]->SetNsigmaTPCle250(2.0);
         dtc1etaphitpc[aniter1]->SetNsigmaTPC250_400(2.0);
          dtc1etaphitpc[aniter1]->SetNsigmaTPC400_450(2.0);
          dtc1etaphitpc[aniter1]->SetNsigmaTPC450_500(1.5); 
          dtc1etaphitpc[aniter1]->SetNsigmaTOF500_800(2.0);
          dtc1etaphitpc[aniter1]->SetNsigmaTOF800_1000(1.5);
          dtc1etaphitpc[aniter1]->SetNsigmaTOFge1000(1.0);
          dtc1etaphitpc[aniter1]->SetRemoveKinks(kTRUE);
          dtc1etaphitpc[aniter1]->SetLabel(kFALSE);
    //----------------------2particle----------------------< KR 
             dtc2etaphitpc[aniter1]=new AliFemtoKKTrackCut();
            dtc2etaphitpc[aniter1]->SetCharge(-1.0);
            dtc2etaphitpc[aniter1]->SetPt(0.14,0.45);
            dtc2etaphitpc[aniter1]->SetEta(-0.8,0.8);
          //PID method
            dtc2etaphitpc[aniter1]->SetMass(KaonMass);
            dtc2etaphitpc[aniter1]->SetMostProbableKaon();
            
          dtc2etaphitpc[aniter1]->SetNsigmaTPCle250(2.0);
         dtc2etaphitpc[aniter1]->SetNsigmaTPC250_400(2.0);
          dtc2etaphitpc[aniter1]->SetNsigmaTPC400_450(2.0);
          dtc2etaphitpc[aniter1]->SetNsigmaTPC450_500(1.5);
          dtc2etaphitpc[aniter1]->SetNsigmaTPCge500(3.0);
          dtc2etaphitpc[aniter1]->SetNsigmaTOF500_800(2.0);
          dtc2etaphitpc[aniter1]->SetNsigmaTOF800_1000(1.5);
          dtc2etaphitpc[aniter1]->SetNsigmaTOFge1000(1.0);
                dtc2etaphitpc[aniter1]->SetRemoveKinks(kTRUE);
//end of strong cuts--------------------------------------------------------------------          	  
          dtc2etaphitpc[aniter1]->SetLabel(kFALSE);
          cutPass1YPtetaphitpc[aniter1] = new AliFemtoCutMonitorParticleYPt(Form("cut2Pass1%stpcM%i", chrgs[ichg], imult), 0.493677);
          cutFail1YPtetaphitpc[aniter1] = new AliFemtoCutMonitorParticleYPt(Form("cut2Fail1%stpcM%i", chrgs[ichg], imult), 0.493677);
            dtc1etaphitpc[aniter1]->AddCutMonitor(cutPass1YPtetaphitpc[aniter1], cutFail1YPtetaphitpc[aniter1]);
            
            cutPass2YPtetaphitpc[aniter1] = new AliFemtoCutMonitorParticleYPt(Form("cut2Pass2%stpcM%i", chrgs[ichg+1], imult), 0.493677); //ichg+1 --> ichg=1 --> charge = -1
            cutFail2YPtetaphitpc[aniter1] = new AliFemtoCutMonitorParticleYPt(Form("cut2Fail2%stpcM%i", chrgs[ichg+1], imult), 0.493677);
          
          dtc2etaphitpc[aniter1]->AddCutMonitor(cutPass2YPtetaphitpc[aniter1], cutFail2YPtetaphitpc[aniter1]);
/*****************************************************/

          cutPass1PIDetaphitpc[aniter1] = new AliFemtoCutMonitorParticlePID(Form("cut2Pass1%stpcM%i", chrgs[ichg], imult),1);
          cutFail1PIDetaphitpc[aniter1] = new AliFemtoCutMonitorParticlePID(Form("cut2Fail1%stpcM%i", chrgs[ichg], imult),1);
          dtc1etaphitpc[aniter1]->AddCutMonitor(cutPass1PIDetaphitpc[aniter1], cutFail1PIDetaphitpc[aniter1]);
 
          sqpcetaphitpc[aniter1] = new AliFemtoPairCutAntiGamma();
       
/*****************************************************/
          anetaphitpc1[aniter1]->SetEventCut(mecetaphitpc[aniter1]);
          anetaphitpc1[aniter1]->SetFirstParticleCut(dtc1etaphitpc[aniter1]);
          anetaphitpc1[aniter1]->SetSecondParticleCut(dtc2etaphitpc[aniter1]); //druga czastka
          anetaphitpc1[aniter1]->SetPairCut(sqpcetaphitpc[aniter1]); 
/*****************************************************/
          //Qinv (without kT bins)
          cqinvkttpc[aniter1] = new AliFemtoQinvCorrFctn(Form("2cqinv%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
          anetaphitpc1[aniter1]->AddCorrFctn(cqinvkttpc[aniter1]);

          //3D cartesian (without kT bins)
          if(run3d){
            cq3dlcmskttpc[aniter1] = new AliFemtoCorrFctn3DLCMSSym(Form("2cq3d%stpcM%i", chrgs[ichg], imult),100,0.5);
            anetaphitpc1[aniter1]->AddCorrFctn(cq3dlcmskttpc[aniter1]);
          }

/*****************************************************/
          if (runktdep) {
            int ktm;
            for (int ikt=0; ikt<2; ikt++) {
              ktm = aniter1*2 + ikt;
              ktpcuts[ktm] = new AliFemtoKTPairCut(ktrng[ikt], ktrng[ikt+1]);
/*****************************************************/
              cqinvkttpc[ktm] = new AliFemtoQinvCorrFctn(Form("2cqinv%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0, shqmax);
              cqinvkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
              anetaphitpc1[aniter1]->AddCorrFctn(cqinvkttpc[ktm]);

              cqinvsqtpc[ktm] = new AliFemtoShareQualityCorrFctn(Form("2cqinvsq%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
              cqinvsqtpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
              anetaphitpc1[aniter1]->AddCorrFctn(cqinvsqtpc[ktm]);
/*****************************************************/
              cqinvinnertpc[ktm] = new AliFemtoTPCInnerCorrFctn(Form("2cqinvinner%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
              cqinvinnertpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
              cqinvinnertpc[ktm]->SetRadius(1.6);
              anetaphitpc1[aniter1]->AddCorrFctn(cqinvinnertpc[ktm]);
              cgamma[aniter1] = new AliFemtoCorrFctnGammaMonitor(Form("2cgammaM%ikT%i", imult, ikt),200,200);
              anetaphitpc1[aniter1]->AddCorrFctn(cgamma[aniter1]);

              if (run3d) {
            cq3dlcmskttpc[ktm] = new AliFemtoCorrFctn3DLCMSSym(Form("2cq3d%stpcM%ikT%i", chrgs[ichg], imult, ikt),60,0.5);
            cq3dlcmskttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
            anetaphitpc1[aniter]->AddCorrFctn(cq3dlcmskttpc[ktm]);
              }
            }
          }

          cdedpetaphi[aniter1] = new AliFemtoCorrFctnDEtaDPhi(Form("2cdedp%stpcM%i", chrgs[ichg], imult),39, 39);
          anetaphitpc1[aniter1]->AddCorrFctn(cdedpetaphi[aniter1]);

          Manager->AddAnalysis(anetaphitpc1[aniter1]);	
        }
    }
  }
//===========================================================end of strong cuts analysis============================================
//===========================================================analysis withouy gamma cut=================================================
  for (int imult=0; imult<3; imult++) {
    if (runmults[imult]) {
        if (runch[ichg]) {
          aniter2 = ichg*3+imult;

          anetaphitpc2[aniter2] = new AliFemtoVertexMultAnalysis(10, -10.0, 10.0, 4, multbins[imult], multbins[imult+1]);
          anetaphitpc2[aniter2]->SetNumEventsToMix(5);
          anetaphitpc2[aniter2]->SetMinSizePartCollection(1);
          anetaphitpc2[aniter2]->SetVerboseMode(kTRUE); 

          mecetaphitpc[aniter2] = new AliFemtoBasicEventCut();
          mecetaphitpc[aniter2]->SetEventMult(0,10000);
          mecetaphitpc[aniter2]->SetVertZPos(-10,10);

          cutPassEvMetaphitpc[aniter2] = new AliFemtoCutMonitorEventMult(Form("cut3Pass%stpcM%i", chrgs[ichg], imult));
          cutFailEvMetaphitpc[aniter2] = new AliFemtoCutMonitorEventMult(Form("cut3Fail%stpcM%i", chrgs[ichg], imult));
          mecetaphitpc[aniter2]->AddCutMonitor(cutPassEvMetaphitpc[aniter2], cutFailEvMetaphitpc[aniter2]);

          cutPassEvVetaphitpc[aniter2] = new AliFemtoCutMonitorEventVertex(Form("cut3Pass%stpcM%i", chrgs[ichg], imult));
          cutFailEvVetaphitpc[aniter2] = new AliFemtoCutMonitorEventVertex(Form("cut3Fail%stpcM%i", chrgs[ichg], imult));
          mecetaphitpc[aniter2]->AddCutMonitor(cutPassEvVetaphitpc[aniter2], cutFailEvVetaphitpc[aniter2]);

          cutPassColletaphitpc[aniter2] = new AliFemtoCutMonitorCollections(Form("cut3Pass%stpcM%i", chrgs[ichg], imult));
              cutFailColletaphitpc[aniter2] = new AliFemtoCutMonitorCollections(Form("cut3Fail%stpcM%i", chrgs[ichg], imult));
              mecetaphitpc[aniter2]->AddCutMonitor(cutPassColletaphitpc[aniter2], cutFailColletaphitpc[aniter2]);
            
            
//-----------------------1 particle-------------------------------------------<
            dtc1etaphitpc[aniter2] = new AliFemtoKKTrackCut();
            dtc1etaphitpc[aniter2]->SetCharge(1.0);
            dtc1etaphitpc[aniter2]->SetPt(0.14,1.5);
            dtc1etaphitpc[aniter2]->SetEta(-0.8,0.8);
          //PID method
            dtc1etaphitpc[aniter2]->SetMass(KaonMass);
            dtc1etaphitpc[aniter2]->SetMostProbableKaon();

         dtc1etaphitpc[aniter2]->SetNsigmaTPCle250(2.0);
         dtc1etaphitpc[aniter2]->SetNsigmaTPC250_400(2.0);
          dtc1etaphitpc[aniter2]->SetNsigmaTPC400_450(2.0);
          dtc1etaphitpc[aniter2]->SetNsigmaTPC450_500(2.0);
          dtc1etaphitpc[aniter2]->SetNsigmaTPCge500(3.0); 
          dtc1etaphitpc[aniter2]->SetNsigmaTOF500_800(2.0);
          dtc1etaphitpc[aniter2]->SetNsigmaTOF800_1000(1.5);
          dtc1etaphitpc[aniter2]->SetNsigmaTOFge1000(1.0);
          //------------------- November 2013 ----------------------------------->
        
          dtc1etaphitpc[aniter2]->SetRemoveKinks(kTRUE);
  
          dtc1etaphitpc[aniter2]->SetLabel(kFALSE);
    //----------------------2particle----------------------< KR 
            // dtc2etaphitpc[aniter] = new AliFemtoESDTrackCut();
            dtc2etaphitpc[aniter2]=new AliFemtoKKTrackCut();
            dtc2etaphitpc[aniter2]->SetCharge(-1.0);
            dtc2etaphitpc[aniter2]->SetPt(0.14,1.5);
            dtc2etaphitpc[aniter2]->SetEta(-0.8,0.8);
          //PID method
            dtc2etaphitpc[aniter2]->SetMass(KaonMass);
            dtc2etaphitpc[aniter2]->SetMostProbableKaon();
          //dtc2etaphitpc[aniter]->SetPIDMethod(AliFemtoESDTrackCut::kContour);
            //------------------- November 2013 -----------------------------------< 
          // new cuts to remove electron (do not take into analysis if 400<p<500) 
         dtc2etaphitpc[aniter2]->SetNsigmaTPCle250(2.0);
         dtc2etaphitpc[aniter2]->SetNsigmaTPC250_400(2.0);
          dtc2etaphitpc[aniter2]->SetNsigmaTPC400_450(2.0);
          dtc2etaphitpc[aniter2]->SetNsigmaTPC450_500(2.0);
          dtc2etaphitpc[aniter2]->SetNsigmaTPCge500(3.0);    
          // new cuts are stronger, better separation of pion in TOF 
          // when momentum is greater then 800 MeV/c
          dtc2etaphitpc[aniter2]->SetNsigmaTOF500_800(2.0);
          dtc2etaphitpc[aniter2]->SetNsigmaTOF800_1000(1.5);
          dtc2etaphitpc[aniter2]->SetNsigmaTOFge1000(1.0);
          //------------------- November 2013 ----------------------------------->
          //Track quality cuts
      
          dtc2etaphitpc[aniter2]->SetRemoveKinks(kTRUE);
          	  
          dtc2etaphitpc[aniter2]->SetLabel(kFALSE);
          cutPass1YPtetaphitpc[aniter2] = new AliFemtoCutMonitorParticleYPt(Form("cut3Pass1%stpcM%i", chrgs[ichg], imult), 0.493677);
          cutFail1YPtetaphitpc[aniter2] = new AliFemtoCutMonitorParticleYPt(Form("cut3Fail1%stpcM%i", chrgs[ichg], imult), 0.493677);
            dtc1etaphitpc[aniter2]->AddCutMonitor(cutPass1YPtetaphitpc[aniter2], cutFail1YPtetaphitpc[aniter2]);
            
            cutPass2YPtetaphitpc[aniter2] = new AliFemtoCutMonitorParticleYPt(Form("cut3Pass2%stpcM%i", chrgs[ichg+1], imult), 0.493677); //ichg+1 --> ichg=1 --> charge = -1
            cutFail2YPtetaphitpc[aniter2] = new AliFemtoCutMonitorParticleYPt(Form("cut3Fail2%stpcM%i", chrgs[ichg+1], imult), 0.493677);
          
          dtc2etaphitpc[aniter2]->AddCutMonitor(cutPass2YPtetaphitpc[aniter2], cutFail2YPtetaphitpc[aniter2]);
/*****************************************************/

          cutPass1PIDetaphitpc[aniter2] = new AliFemtoCutMonitorParticlePID(Form("cut3Pass1%stpcM%i", chrgs[ichg], imult),1);
          cutFail1PIDetaphitpc[aniter2] = new AliFemtoCutMonitorParticlePID(Form("cut3Fail1%stpcM%i", chrgs[ichg], imult),1);
          dtc1etaphitpc[aniter2]->AddCutMonitor(cutPass1PIDetaphitpc[aniter2], cutFail1PIDetaphitpc[aniter2]);
 
          sqpcetaphitpc[aniter2] = new AliFemtoPairCutAntiGamma();
       
/*****************************************************/
          anetaphitpc2[aniter2]->SetEventCut(mecetaphitpc[aniter2]);
          anetaphitpc2[aniter2]->SetFirstParticleCut(dtc1etaphitpc[aniter2]);
          anetaphitpc2[aniter2]->SetSecondParticleCut(dtc2etaphitpc[aniter2]); //druga czastka
          anetaphitpc2[aniter2]->SetPairCut(sqpcetaphitpc[aniter2]); 
/*****************************************************/
          //Qinv (without kT bins)
          cqinvkttpc[aniter2] = new AliFemtoQinvCorrFctn(Form("3cqinv%stpcM%i", chrgs[ichg], imult),nbinssh,0.0,shqmax);
          anetaphitpc2[aniter2]->AddCorrFctn(cqinvkttpc[aniter2]);

          //3D cartesian (without kT bins)
          if(run3d){
            cq3dlcmskttpc[aniter2] = new AliFemtoCorrFctn3DLCMSSym(Form("3cq3d%stpcM%i", chrgs[ichg], imult),100,0.5);
            anetaphitpc2[aniter2]->AddCorrFctn(cq3dlcmskttpc[aniter2]);
          }

/*****************************************************/
          if (runktdep) {
            int ktm;
            for (int ikt=0; ikt<2; ikt++) {
              ktm = aniter2*2 + ikt;
              ktpcuts[ktm] = new AliFemtoKTPairCut(ktrng[ikt], ktrng[ikt+1]);
/*****************************************************/
              cqinvkttpc[ktm] = new AliFemtoQinvCorrFctn(Form("3cqinv%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0, shqmax);
              cqinvkttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
              anetaphitpc2[aniter2]->AddCorrFctn(cqinvkttpc[ktm]);

              cqinvsqtpc[ktm] = new AliFemtoShareQualityCorrFctn(Form("3cqinvsq%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
              cqinvsqtpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
              anetaphitpc2[aniter2]->AddCorrFctn(cqinvsqtpc[ktm]);
/*****************************************************/
              cqinvinnertpc[ktm] = new AliFemtoTPCInnerCorrFctn(Form("3cqinvinner%stpcM%ikT%i", chrgs[ichg], imult, ikt),nbinssh,0.0,shqmax);
              cqinvinnertpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
              cqinvinnertpc[ktm]->SetRadius(1.6);
              anetaphitpc2[aniter2]->AddCorrFctn(cqinvinnertpc[ktm]);
              cgamma[aniter2] = new AliFemtoCorrFctnGammaMonitor(Form("3cgammaM%ikT%i", imult, ikt),200,200);
              anetaphitpc2[aniter2]->AddCorrFctn(cgamma[aniter2]);

              if (run3d) {
            cq3dlcmskttpc[ktm] = new AliFemtoCorrFctn3DLCMSSym(Form("3cq3d%stpcM%ikT%i", chrgs[ichg], imult, ikt),60,0.5);
            cq3dlcmskttpc[ktm]->SetPairSelectionCut(ktpcuts[ktm]);
            anetaphitpc2[aniter2]->AddCorrFctn(cq3dlcmskttpc[ktm]);
              }
            }
          }

          cdedpetaphi[aniter2] = new AliFemtoCorrFctnDEtaDPhi(Form("3cdedp%stpcM%i", chrgs[ichg], imult),39, 39);
          anetaphitpc2[aniter2]->AddCorrFctn(cdedpetaphi[aniter2]);

          Manager->AddAnalysis(anetaphitpc2[aniter2]);	
        }
    }
  }
 //=============================================================end of analysis without gamma cut====================================================

  return Manager;
}                         
                      
