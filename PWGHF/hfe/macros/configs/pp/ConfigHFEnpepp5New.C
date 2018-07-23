TF1* GetTPCCorrection(TString type){
   TString Map="TPCCorrMaps_5TeV.root";
   TFile *f = TFile::Open(Form("$ALICE_PHYSICS/PWGHF/hfe/macros/%s", Map.Data()));
   if (!f->IsOpen()) {
      printf("TPC map file not found \n");
      return 0;
   }
   gROOT->cd();
   TF1* CorrectionFunction(NULL);
   if (type == "Momentum") {
      CorrectionFunction = (TF1*) f->GetObjectChecked("TPCnsigmavspcor","TF1");
      if(!CorrectionFunction) printf("TPCCorrMaps_5TeV.root:TPCnsigmavspcor not found \n");
   } else if (type == "Momentum_Width"){
      CorrectionFunction = (TF1*) f->GetObjectChecked("TPCnsigmavspcor_width","TF1");
      if(!CorrectionFunction) printf("TPCCorrMaps_5TeV.root:TPCnsigmavspcor_width not found \n");
   } else if (type == "Eta"){
      CorrectionFunction = (TF1*) f->GetObjectChecked("TPCnsigmavsetacor","TF1");
      if(!CorrectionFunction) printf("TPCCorrMaps_5TeV.root:TPCnsigmavsetacor not found \n");
   } else if (type == "Eta_Width"){
      CorrectionFunction = (TF1*) f->GetObjectChecked("TPCnsigmavsetacor_width","TF1");
      if(!CorrectionFunction) printf("TPCCorrMaps_5TeV.root:TPCnsigmavsetacor_width not found \n");
   }
   return CorrectionFunction;
}

AliAnalysisTaskHFE* ConfigHFEnpepp5New(Bool_t useMC, Bool_t isAOD, TString appendix,
                                    UChar_t TPCcl=70, UChar_t TPCclPID = 80,
                                    UChar_t ITScl=3, Double_t DCAxy=1000., Double_t DCAz=1000.,
                                    Double_t* tpcdEdxcutlow=NULL, Double_t* tpcdEdxcuthigh=NULL,
                                    Double_t TOFs=3., Int_t TOFmis=0,
                                    Int_t itshitpixel = 0, Int_t ikink = 0,
                                    Double_t etami=-0.8, Double_t etama=0.8,
                                    Double_t assETAm=-0.8, Double_t assETAp=0.8,
                                    Double_t assMinPt=0.2, Int_t assITS=2,
                                    Int_t assTPCcl=100, Int_t assTPCPIDcl=80,
                                    Double_t assDCAr=1.0, Double_t assDCAz=2.0,
                                    Double_t *assTPCSminus=NULL, Double_t *assTPCSplus=NULL,
                                    Double_t assITSpid=-3., Double_t assTOFs=3.,
                                    Bool_t useCat1Tracks = kTRUE, Bool_t useCat2Tracks = kTRUE, Int_t weightlevelback = -1,
                                    Int_t HadronContFunc=0, Int_t PrimaryVertexTyp)
{

   Bool_t etacor = false;
   Bool_t MomCor = false;

   //***************************************//
   //        Setting up the HFE cuts        //
   //***************************************//

   AliHFEcuts *hfecuts = new AliHFEcuts(appendix,"HFE cuts for pp");
   //hfecuts->SetQAOn();
   hfecuts->CreateStandardCuts();
   if(isAOD) hfecuts->SetAODFilterBit(4); // standard cuts with very loose DCA

   hfecuts->SetMinNClustersTPC(TPCcl);
   hfecuts->SetMinNClustersTPCPID(TPCclPID);
   hfecuts->SetMinNClustersITS(ITScl);
   hfecuts->SetMinRatioTPCclusters(0.6);
   hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
   hfecuts->SetCutITSpixel(itshitpixel);
   hfecuts->SetCheckITSLayerStatus(kFALSE);
   hfecuts->SetEtaRange(etami,etama);
   hfecuts->SetRejectKinkDaughters();
   if (ikink == 0){
      hfecuts->SetRejectKinkMothers();
   } else {
      hfecuts->SetAcceptKinkMothers();
   }
   hfecuts->SetMaxImpactParam(DCAxy,DCAz);

   if(PrimaryVertexTyp==1){
      hfecuts->SetUseMixedVertex(kFALSE);
      hfecuts->SetUseSPDVertex(kFALSE);
      hfecuts->SetUseTrackVertex(kTRUE);
   } else if(PrimaryVertexTyp==2){
      hfecuts->SetUseMixedVertex(kFALSE);
      hfecuts->SetUseSPDVertex(kTRUE);
      hfecuts->SetUseTrackVertex(kFALSE);
   } else{
      hfecuts->SetUseMixedVertex(kTRUE);
      hfecuts->SetUseSPDVertex(kFALSE);
      hfecuts->SetUseTrackVertex(kFALSE);
   }

   hfecuts->SetVertexRange(10.);

   if (hfecuts->GetUseSPDVertex() || hfecuts->GetUseMixedVertex() ) {
      //if SPD vertices are used acivate resolution cut
      printf("########## Additional cut for SPD vertex resolution ######### \n");
      hfecuts->SetSPDVtxResolutionCut();
   }


   // TOF settings:
   Int_t usetof=0;
   Bool_t kTOFmis=kFALSE;
   if (TOFs>0.){
      usetof = 1;
      printf("CONFIGURATION FILE: TOF is used \n");
      hfecuts->SetTOFPIDStep(kTRUE);
      printf("CONFIGURATION FILE: TOF PID step is requested !!!! \n");
      if (TOFmis>0){
         kTOFmis = kTRUE;
         printf("CONFIGURATION FILE: TOF mismatch rejection is set ON \n");
      }
   }

   //***************************************//
   //        Setting up the task            //
   //***************************************//

   AliAnalysisTaskHFE *task = new AliAnalysisTaskHFE(Form("HFEtask%s",appendix.Data()));
   printf("task %p\n", task);
   task->SetppAnalysis();
   //if(!isAOD) task->SetRemoveFirstEventInChunk();
   task->SetRemovePileUp(kTRUE);
   task->SetHFECuts(hfecuts);
   task->GetPIDQAManager()->SetHighResolutionHistos();

   //***************************************//
   //          Variable manager             //
   //***************************************//
   // Define Variables
   //standard binning
   Double_t ptbinning[36] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10., 12., 14., 16., 18., 20.};
   Double_t etabinning[17] = {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

   Double_t phibinning[5] = {0, 0.5*TMath::Pi(), TMath::Pi(), 1.5*TMath::Pi(),2*TMath::Pi()};

   Int_t sizept=(sizeof(ptbinning)/sizeof(double))-1;
   Int_t sizeeta=(sizeof(etabinning)/sizeof(double))-1;
   Int_t sizephi=(sizeof(phibinning)/sizeof(double))-1;

   AliHFEvarManager *vm = task->GetVarManager();
   vm->AddVariable("pt", sizept, ptbinning);
   vm->AddVariable("eta", sizeeta, -0.8,0.8);
   vm->AddVariable("phi", sizephi, phibinning);
   vm->AddVariable("charge");
   vm->AddVariable("source");

   // The part dedicated to the background subtraction
   // should be implemented in a different way, reading it from a root file.

   if(!useMC){

      enum HadronContaminationFunctionChoice
      {
         kHadronDown = 1,
         kHadronUp= 2,
         kHadronError = 3,

         kHadronTPCMinus05 = 4,
         kHadronTPC0 = 5,
         kHadronTPC025 = 6,
         kHadronTPCMinus025 = 7,
         kHadronTPCMinus075 = 8,
         kHadronTPCMinus13 = 9,
         kHadronTPCUp10 = 10,
         kHadronTPCUp15 = 11,
         kHadronTPCUp20 = 12,
         kHadronEtaRange07 = 13,
         kHadronEtaRange06 = 14,
         kHadronEtaRange05 = 15,
         kHadronTPCMinus1_Old = 16,
         
         kHadronTPConly_Minus05 = 20,
         kHadronTPConly_Minus025 = 21,
         kHadronTPConly_Minus0125 = 22,
         kHadronTPConly_0125 = 23,
         kHadronTPConly_025 = 24,
         kHadronTPConly_05 = 25,
         kHadronTPConly_Landau = 26,
         kHadronTPConly_Down = 27,
         kHadronTPConly_Up = 28,
         kHadronTPConly_EtaRange07 = 29,
         kHadronTPConly_EtaRange06 = 30,
         kHadronTPConly_EtaRange05 = 31,
         kHadronTPConly_0_Old = 32


         
      };

      // First hadron contamination fit for 5TeV  by Sebastian Hornung, March 9, 2017
      // relative to the case of a TPC PID cut at -1 sigma

      TF1 *hBackground;
      if(TOFs==0){
         // Default of the TPC-TOF - TO-BE-UPDATED
         switch (HadronContFunc) {
                case kHadronTPConly_Minus05:
               TF1 *hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])+ [3]*TMath::Gaus(x,[4],[5])",0.,60.);
               hBackground->SetParameter(0,5.14946e+00);
               hBackground->SetParameter(1,1.58931e+01);
               hBackground->SetParameter(2,4.03552e+00);

               hBackground->SetParameter(3,1.70000e-01);
               hBackground->SetParameter(4,1.85645e+00);
               hBackground->SetParameter(5,1.20130e-01);

               break;
               
                case kHadronTPConly_Minus025:
               TF1 *hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])+ [3]*TMath::Gaus(x,[4],[5])",0.,60.);
               hBackground->SetParameter(0,4.29739e+00);
               hBackground->SetParameter(1,1.80873e+01);
               hBackground->SetParameter(2,4.61691e+00);

               hBackground->SetParameter(3,1.70000e-01);
               hBackground->SetParameter(4,1.84226e+00);
               hBackground->SetParameter(5,1.18907e-01);

               break;
               
                case kHadronTPConly_Minus0125:
               TF1 *hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])+ [3]*TMath::Gaus(x,[4],[5])",0.,60.);
               hBackground->SetParameter(0,3.28834e+00);
               hBackground->SetParameter(1,1.86554e+01);
               hBackground->SetParameter(2,4.76731e+00);

               hBackground->SetParameter(3,1.70000e-01);
               hBackground->SetParameter(4,1.83479e+00);
               hBackground->SetParameter(5,1.18766e-01);

               break;
               
                case kHadronTPConly_0125:
                    
               TF1 *hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])+ [3]*TMath::Gaus(x,[4],[5])",0.,60.);
               hBackground->SetParameter(0,1.76294e+00);
               hBackground->SetParameter(1,1.96221e+01);
               hBackground->SetParameter(2,5.03374e+00);

               hBackground->SetParameter(3,1.70000e-01);
               hBackground->SetParameter(4,1.81881e+00);
               hBackground->SetParameter(5,1.18879e-01);


               break;
               
                case kHadronTPConly_025:
                    
               TF1 *hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])+ [3]*TMath::Gaus(x,[4],[5])",0.,60.);
               hBackground->SetParameter(0,1.13397e+00);
               hBackground->SetParameter(1,1.96791e+01);
               hBackground->SetParameter(2,5.06089e+00);

               hBackground->SetParameter(3,1.70000e-01);
               hBackground->SetParameter(4,1.81024e+00);
               hBackground->SetParameter(5,1.19064e-01);


                    
               break;
               
                case kHadronTPConly_05:
                    
               hBackground = new TF1("hadronicBackgroundFunction", "[0]+[1]*TMath::Erf([2]*x+[3]) + [4]*TMath::Gaus(x,[5],[6])",0.,60.);
               hBackground->SetParameter(0,3.16386e-01);
               hBackground->SetParameter(1,3.16432e-01);
               hBackground->SetParameter(2,1.96531e-01);
               hBackground->SetParameter(3,-3.15578e+00);
               hBackground->SetParameter(4,1.70000e-01);
               hBackground->SetParameter(5,1.79193e+00);
               hBackground->SetParameter(6,1.19514e-01);

                    
               break;
               
               case kHadronTPConly_Landau:

               hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])",0. ,60.);
               hBackground->SetParameter(0,2.00832e+00);
               hBackground->SetParameter(1,1.82135e+01);
               hBackground->SetParameter(2,4.62951e+00);

               break;

                case kHadronTPConly_Down:

               hBackground = new TF1("hadronicBackgroundFunction", "[0]+[1]*TMath::Erf([2]*x+[3])",0. ,60.);
               hBackground->SetParameter(0,2.56382e-01);
               hBackground->SetParameter(1,2.56453e-01);
               hBackground->SetParameter(2,3.10993e-01);
               hBackground->SetParameter(3,-3.15362e+00);

               break;

                case kHadronTPConly_Up:

               hBackground = new TF1("hadronicBackgroundFunction", "[0]+[1]*TMath::Erf([2]*x+[3])",0. ,60.);
               hBackground->SetParameter(0,2.14571e-02);
               hBackground->SetParameter(1,2.14768e-02);
               hBackground->SetParameter(2,4.87058e-01);
               hBackground->SetParameter(3,-4.03708e+00);

               break;
               
                case kHadronTPConly_EtaRange07:

               hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])+ [3] * TMath::Gaus(x, [4], [5])", 0., 30);
               hBackground->SetParameter(0,3.75777e+00);
               hBackground->SetParameter(1,2.03415e+01);
               hBackground->SetParameter(2,5.20870e+00);

               hBackground->SetParameter(3,1.70000e-01);
               hBackground->SetParameter(4,1.84239e+00);
               hBackground->SetParameter(5,1.28382e-01);

               break;
               
                case kHadronTPConly_EtaRange06:

               hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])+ [3] * TMath::Gaus(x, [4], [5])", 0., 30);
               hBackground->SetParameter(0,3.02246e+00);
               hBackground->SetParameter(1,1.85956e+01);
               hBackground->SetParameter(2,4.69020e+00);

               hBackground->SetParameter(3,1.70000e-01);
               hBackground->SetParameter(4,1.81410e+00);
               hBackground->SetParameter(5,1.49956e-01);


               break;
               
                case kHadronTPConly_EtaRange05:

               hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])+ [3] * TMath::Gaus(x, [4], [5])", 0., 30);
               hBackground->SetParameter(0,2.50954e+00);
               hBackground->SetParameter(1,1.70105e+01);
               hBackground->SetParameter(2,4.21313e+00);
 
               hBackground->SetParameter(3,1.70000e-01);
               hBackground->SetParameter(4,1.81304e+00);
               hBackground->SetParameter(5,1.74524e-01);


               break;

                case kHadronTPConly_0_Old:
               hBackground = new TF1("hadronicBackgroundFunction", "[0]+[1]*TMath::Erf([2]*x+[3])",0. ,30.);

               hBackground->SetParameter(0,1.29889e-01);
               hBackground->SetParameter(1,1.29923e-01);
               hBackground->SetParameter(2,3.31964e-01);
               hBackground->SetParameter(3,-3.23896e+00);

               break;



     
                default:
                    
               hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])+ [3] * TMath::Gaus(x, [4], [5])", 0., 30);
               hBackground->SetParameter(0,2.55254e+00);
               hBackground->SetParameter(1,1.93266e+01);
               hBackground->SetParameter(2,4.94811e+00);

               hBackground->SetParameter(3,1.70000e-01);
               hBackground->SetParameter(4,1.82700e+00);
               hBackground->SetParameter(5,1.18751e-01);


               break;
               
               


         }
      }else{
         switch (HadronContFunc) {
            case kHadronDown:
               hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])+ [3] * TMath::Gaus(x, [4], [5])", 0., 30);
               hBackground->SetParameter(0,2.40611e+00);
               hBackground->SetParameter(1,1.23231e+01);
               hBackground->SetParameter(2,3.11748e+00);

               hBackground->SetParameter(3,2.53125e-02);
               hBackground->SetParameter(4,8.72447e-01);
               hBackground->SetParameter(5,2.94525e-02);
               break;
            case kHadronUp:
               hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])+ [3] * TMath::Gaus(x, [4], [5])", 0., 30);
               hBackground->SetParameter(0,5.30204e-01);
               hBackground->SetParameter(1,9.28479e+00);
               hBackground->SetParameter(2,2.11948e+00);

               hBackground->SetParameter(3,2.51495e-02);
               hBackground->SetParameter(4,8.72771e-01);
               hBackground->SetParameter(5,2.82148e-02);
               break;
            case kHadronError:
               hBackground = new TF1("hadronicBackgroundFunction", "[0]+[1]*TMath::Erf([2]*x+[3]) + [4] * TMath::Gaus(x, [5], [6])",0. ,60.);
               hBackground->SetParameter(0,4.98687e-01);
               hBackground->SetParameter(1,4.98692e-01);
               hBackground->SetParameter(2,5.46485e-01);
               hBackground->SetParameter(3,-3.45804e+00);

               hBackground->SetParameter(4,2.00000e-02);
               hBackground->SetParameter(5,8.78696e-01);
               hBackground->SetParameter(6,2.62630e-02);
               break;
            case kHadronTPC025:
               hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])+ [3] * TMath::Gaus(x, [4], [5])", 0., 30);
               hBackground->SetParameter(0,4.27361e+00);
               hBackground->SetParameter(1,1.90821e+01);
               hBackground->SetParameter(2,4.75257e+00);

               //remaining Kaon crossing
               hBackground->SetParameter(3,3.00000e-02);
               hBackground->SetParameter(4,8.78363e-01);
               hBackground->SetParameter(5,2.72400e-02);
               break;
            case kHadronTPC0:
               hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])+ [3] * TMath::Gaus(x, [4], [5])", 0., 30);
               hBackground->SetParameter(0,4.84277e+00);
               hBackground->SetParameter(1,1.72791e+01);
               hBackground->SetParameter(2,4.28400e+00);

               //remaining Kaon crossing
               hBackground->SetParameter(3,3.00000e-02);
               hBackground->SetParameter(4,8.77877e-01);
               hBackground->SetParameter(5,2.62314e-02);
               break;
            case kHadronTPCMinus025:
               hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])+ [3] * TMath::Gaus(x, [4], [5])", 0., 30);
               hBackground->SetParameter(0,4.93950e+00);
               hBackground->SetParameter(1,1.53102e+01);
               hBackground->SetParameter(2,3.77895e+00);

               //remaining Kaon crossing
               hBackground->SetParameter(3,2.00000e-02);
               hBackground->SetParameter(4,8.79282e-01);
               hBackground->SetParameter(5,2.79708e-02);
               break;
            case kHadronTPCMinus075:
               hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])+ [3] * TMath::Gaus(x, [4], [5])", 0., 30);
               hBackground->SetParameter(0,6.20530e+00);
               hBackground->SetParameter(1,1.25144e+01);
               hBackground->SetParameter(2,3.09613e+00);

               //remaining Kaon crossing
               hBackground->SetParameter(3,2.00000e-02);
               hBackground->SetParameter(4,8.78686e-01);
               hBackground->SetParameter(5,2.67958e-02);
               break;
            case kHadronTPCMinus05:
               hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])+ [3] * TMath::Gaus(x, [4], [5])", 0., 30);
               hBackground->SetParameter(0,6.62688e+00);
               hBackground->SetParameter(1,1.42852e+01);
               hBackground->SetParameter(2,3.53148e+00);

               //remaining Kaon crossing
               hBackground->SetParameter(3,2.00000e-02);
               hBackground->SetParameter(4,8.78893e-01);
               hBackground->SetParameter(5,2.73265e-02);
               break;
            case kHadronTPCMinus13:
               hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])+ [3] * TMath::Gaus(x, [4], [5])", 0., 30);
               hBackground->SetParameter(0,5.21317e+00);
               hBackground->SetParameter(1,9.14025e+00);
               hBackground->SetParameter(2,2.25821e+00);

               //remaining Kaon crossing
               hBackground->SetParameter(3,2.00000e-02);
               hBackground->SetParameter(4,8.78584e-01);
               hBackground->SetParameter(5,-2.60403e-02);
               break;
            case kHadronTPCUp10:
               hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])+ [3] * TMath::Gaus(x, [4], [5])", 0., 30);
               hBackground->SetParameter(0,5.84601e+00);
               hBackground->SetParameter(1,1.05731e+01);
               hBackground->SetParameter(2,2.60824e+00);

               hBackground->SetParameter(3,5.00000e-03);
               hBackground->SetParameter(4,8.84228e-01);
               hBackground->SetParameter(5,2.96319e-02);
               break;
            case kHadronTPCUp15:
               hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])+ [3] * TMath::Gaus(x, [4], [5])", 0., 30);
               hBackground->SetParameter(0,5.83433e+00);
               hBackground->SetParameter(1,1.07907e+01);
               hBackground->SetParameter(2,2.66667e+00);

               hBackground->SetParameter(3,1.00000e-02);
               hBackground->SetParameter(4,8.81435e-01);
               hBackground->SetParameter(5,2.61541e-02);
               break;
            case kHadronTPCUp20:
               hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])+ [3] * TMath::Gaus(x, [4], [5])", 0., 30);
               hBackground->SetParameter(0,5.67295e+00);
               hBackground->SetParameter(1,1.07832e+01);
               hBackground->SetParameter(2,2.66052e+00);

               hBackground->SetParameter(3,1.00000e-02);
               hBackground->SetParameter(4,8.81985e-01);
               hBackground->SetParameter(5,2.73737e-02);
               break;
            case kHadronEtaRange07:
               hBackground = new TF1("hadronicBackgroundFunction", "[0]+[1]*TMath::Erf([2]*x+[3]) + [4] * TMath::Gaus(x, [5], [6])",0. ,60.);
               hBackground->SetParameter(0,4.99690e-01);
               hBackground->SetParameter(1,4.99690e-01);
               hBackground->SetParameter(2,4.92101e-01);
               hBackground->SetParameter(3,-3.67450e+00);

               hBackground->SetParameter(4,3.22580e-02);
               hBackground->SetParameter(5,8.67630e-01);
               hBackground->SetParameter(6,2.78381e-02);
               break;
            case kHadronEtaRange06:
               hBackground = new TF1("hadronicBackgroundFunction", "[0]+[1]*TMath::Erf([2]*x+[3]) + [4] * TMath::Gaus(x, [5], [6])",0. ,60.);
               hBackground->SetParameter(0,4.99690e-01);
               hBackground->SetParameter(1,4.99690e-01);
               hBackground->SetParameter(2,5.11702e-01);
               hBackground->SetParameter(3,-3.71773e+00);

               hBackground->SetParameter(4,4.29179e-02);
               hBackground->SetParameter(5,8.67045e-01);
               hBackground->SetParameter(6,2.65574e-02);
               break;
            case kHadronEtaRange05:
               hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])+ [3] * TMath::Gaus(x, [4], [5])", 0., 30);
               hBackground->SetParameter(0,6.58853e+01);
               hBackground->SetParameter(1,1.73884e+01);
               hBackground->SetParameter(2,4.21580e+00);

               hBackground->SetParameter(3,3.07715e-02);
               hBackground->SetParameter(4,8.67344e-01);
               hBackground->SetParameter(5,2.84939e-02);
               break;
            case kHadronTPCMinus1_Old:
               hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])+ [3] * TMath::Gaus(x, [4], [5])", 0., 30);
               hBackground->SetParameter(0,6.12980e+00);
               hBackground->SetParameter(1,1.03676e+01);
               hBackground->SetParameter(2,2.59010e+00);

               hBackground->SetParameter(3,3.69613e-02);
               hBackground->SetParameter(4,8.87191e-01);
               hBackground->SetParameter(5,1.47384e-02);
               break;
            default: // lower cut -1
               hBackground = new TF1("hadronicBackgroundFunction", "[0]*TMath::Landau(x,[1],[2])+ [3] * TMath::Gaus(x, [4], [5])", 0., 30);
               hBackground->SetParameter(0,5.80726e+00);
               hBackground->SetParameter(1,1.09115e+01);
               hBackground->SetParameter(2,2.70025e+00);

               hBackground->SetParameter(3,2.00000e-02);
               hBackground->SetParameter(4,8.78616e-01);
               hBackground->SetParameter(5,2.63707e-02);
               break;
         }
      }

      //error function
      task->SetBackGroundFactorsFunction(hBackground);

   }
   
   //***************************************//
   //          Configure the PID            //
   //***************************************//

   // Define PID
   AliHFEpid *pid = task->GetPID();
   if(useMC) pid->SetHasMCData(kTRUE);

   if (usetof){
      pid->AddDetector("TOF", 0);
      pid->AddDetector("TPC", 1);
   } else {
      pid->AddDetector("TPC", 0);
   }

   // Configure TPC PID
   // do the identical thing in data and MC
   Double_t paramsTPCdEdxcutlow[12] ={0.0, 0.0, 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
   if(tpcdEdxcutlow) memcpy(paramsTPCdEdxcutlow,tpcdEdxcutlow,sizeof(paramsTPCdEdxcutlow));

   Double_t paramsTPCdEdxcuthigh[12] ={3.0, 3.0, 3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0};
   if(tpcdEdxcuthigh) memcpy(paramsTPCdEdxcuthigh,tpcdEdxcuthigh,sizeof(paramsTPCdEdxcuthigh));

   char *cutmodel;

   if(useMC){ // constant (default) cut for MC
      cutmodel="pol0(0)";
      Double_t params[1];
      params[0]=paramsTPCdEdxcutlow[0];
      pid->ConfigureTPCdefaultCut(cutmodel, params,tpcdEdxcuthigh[0]);
   } else { // correct for mean shift in data
      cutmodel="pol1(0)+[2]*pol1(3)";
      Double_t params[5];
      params[0]=0.061296;
      params[1]=0.036105;
      params[2]=paramsTPCdEdxcutlow[0];
      params[3]=0.898012;
      params[4]=0.037808;
      pid->ConfigureTPCdefaultCut(cutmodel, params,tpcdEdxcuthigh[0]);
   }

   if(!useMC){
      AliHFEpidTPC *tpcpid = pid->GetDetPID(AliHFEpid::kTPCpid);
      if(etacor){
         printf("CONFIGURATION FILE: Eta correction of electrons in the TPC \n");
         // Apply eta correction
         TF1 *etacorrection = GetTPCCorrection("Eta");
         TF1 *etacorrectionWidth = GetTPCCorrection("Eta_Width");
         if(etacorrection) tpcpid->SetEtaCorrections(etacorrection,etacorrectionWidth);
      }
      if(MomCor){
         printf("CONFIGURATION FILE: Momentum correction of electrons in the TPC \n");
         // Apply eta correction
         TF1 *MomCorrection = GetTPCCorrection("Momentum");
         TF1 *MomCorrectionWidth = GetTPCCorrection("Momentum_Width");
         if(MomCorrection) tpcpid->SetMomentumCorrections(MomCorrection, MomCorrectionWidth);
      }
   }

   // Configure TOF PID
   if (usetof){
      pid->ConfigureTOF(TOFs);
      AliHFEpidTOF *tofpid = pid->GetDetPID(AliHFEpid::kTOFpid);
      if (kTOFmis){
         tofpid->SetRejectTOFmismatch();
      }
   }

   // To make different upper TOF cut to see contamination effect
   // The below two lines should be removed after this check
   //AliHFEpidTOF *tofpid = pid->GetDetPID(AliHFEpid::kTOFpid);
   //if(TOFs<3.) tofpid->SetTOFnSigmaBand(-3,TOFs); //only to check the assymmetric tof cut


   //***************************************//
   //       Configure NPE plugin            //
   //***************************************//

   AliHFENonPhotonicElectron *backe = new AliHFENonPhotonicElectron(Form("HFEBackGroundSubtractionPID2%s",appendix.Data()),"Background subtraction");  //appendix
   //Setting the Cuts for the Associated electron-pool
   AliHFEcuts *hfeBackgroundCuts = new AliHFEcuts(Form("HFEBackSub%s",appendix.Data()),"Background sub Cuts");
   hfeBackgroundCuts->SetEtaRange(assETAm,assETAp);
   hfeBackgroundCuts->SetPtRange(assMinPt,20.);

   hfeBackgroundCuts->SetMaxChi2perClusterTPC(4);

   hfeBackgroundCuts->SetMinNClustersITS(assITS);
   hfeBackgroundCuts->SetMinNClustersTPC(assTPCcl);
   hfeBackgroundCuts->SetMinNClustersTPCPID(assTPCPIDcl);
   hfeBackgroundCuts->SetMaxImpactParam(assDCAr,assDCAz);
   if(isAOD) hfeBackgroundCuts->SetAODFilterBit(0); // Standard TPC only tracks
   hfeBackgroundCuts->SetQAOn();			        // QA

   AliHFEpid *pidbackground = backe->GetPIDBackground();
   if(useMC) pidbackground->SetHasMCData(kTRUE);

   if (assTOFs>0.){
      pidbackground->AddDetector("TOF", 0);
      pidbackground->AddDetector("TPC", 1);
   } else {
      pidbackground->AddDetector("TPC", 0);
   }

   Double_t paramsTPCdEdxcutlowAssoc[12] ={-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0};
   if(assTPCSminus) memcpy(paramsTPCdEdxcutlowAssoc,assTPCSminus,sizeof(paramsTPCdEdxcutlowAssoc));

   Double_t paramsTPCdEdxcuthighAssoc[12] ={3.0, 3.0, 3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0};
   if(assTPCSplus) memcpy(paramsTPCdEdxcuthighAssoc,assTPCSplus,sizeof(paramsTPCdEdxcuthighAssoc));

   char *cutmodelAssoc;
   cutmodelAssoc="pol0";
   for(Int_t a=0;a<11;a++){
      // Not necessary anymore, since the pPb case is handled similarly to the pp case
      //   cout << a << " " << paramsTPCdEdxcut[a] << endl;
      Double_t tpcparamlow[1]={paramsTPCdEdxcutlowAssoc[a]};
      Float_t tpcparamhigh=paramsTPCdEdxcuthighAssoc[a];
      pidbackground->ConfigureTPCcentralityCut(a,cutmodelAssoc,tpcparamlow,tpcparamhigh);
   }
   pidbackground->ConfigureTPCdefaultCut(cutmodelAssoc,paramsTPCdEdxcutlowAssoc,paramsTPCdEdxcuthighAssoc[0]); // After introducing the pPb flag, pPb is merged with pp and this line defines the cut
   //backe->GetPIDBackgroundQAManager()->SetHighResolutionHistos();

   if (assTOFs>0.){
      pidbackground->ConfigureTOF(TOFs);
   }

   backe->SetHFEBackgroundCuts(hfeBackgroundCuts);

   // Selection of associated tracks for the pool
   if(useCat1Tracks) backe->SelectCategory1Tracks(kTRUE);
   if(useCat2Tracks){
      backe->SelectCategory2Tracks(kTRUE);
      backe->SetITSMeanShift(-0.5);
      backe->SetITSnSigmaHigh(assITSpid);
      Double_t assITSminus = -1.0 * assITSpid;
      backe->SetITSnSigmaLow(assITSminus);
      //backe->SetminPt(assMinPt);
   }

   // apply opening angle cut to reduce file size
   backe->SetMaxInvMass(0.6);

   Double_t massbinning[32] = {0.0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.16,0.18,0.2,0.24,0.28,0.32,0.36,0.4,0.44,0.48,0.52,0.56,0.6};
   Int_t sizemass=(sizeof(massbinning)/sizeof(double))-1;

   backe->SetPtBinning(sizept, ptbinning);
   backe->SetEtaBinning(sizeeta, etabinning);
   backe->SetInvMassBinning(sizemass, massbinning);
   backe->SetPhiBinning(sizephi,phibinning);

   if(useMC) {
      if((weightlevelback >=0) && (weightlevelback <3)) backe->SetWithWeights(weightlevelback);
   }

   task->SetHFEBackgroundSubtraction(backe);
   task->SetWeightHist();

   // Also May 30: try to add tagged tracks (photon conversions)
   // not used since long ...
   AliHFEcuts *v0trackCuts = new AliHFEcuts("V0trackCuts", "Track Cuts for tagged track Analysis");
   v0trackCuts->CreateStandardCuts();
   v0trackCuts->SetMinNClustersTPC(TPCcl);
   v0trackCuts->SetMinNClustersTPCPID(TPCclPID);
   v0trackCuts->SetMinRatioTPCclusters(0.6);
   v0trackCuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
   v0trackCuts->SetMinNClustersITS(1);
   v0trackCuts->SetCutITSpixel(AliHFEextraCuts::kFirst);
   v0trackCuts->SetCheckITSLayerStatus(kFALSE);
   v0trackCuts->UnsetVertexRequirement();
   //hfecuts->SetSigmaToVertex(10);
   //v0trackCuts->SetTOFPIDStep(kTRUE);
   v0trackCuts->SetQAOn();

   task->SwitchOnPlugin(AliAnalysisTaskHFE::kTaggedTrackAnalysis);
   task->SetTaggedTrackCuts(v0trackCuts);
   task->SetCleanTaggedTrack(kTRUE);
   // end tagged tracks

   // QA
   printf("task %p\n", task);
   task->SetQAOn(AliAnalysisTaskHFE::kPIDqa);
   task->SetQAOn(AliAnalysisTaskHFE::kMCqa);
   task->SwitchOnPlugin(AliAnalysisTaskHFE::kDEstep);
   task->SwitchOnPlugin(AliAnalysisTaskHFE::kNonPhotonicElectron);

   printf("*************************************\n");
   printf("Configuring standard Task:\n");
   task->PrintStatus();
   pid->PrintStatus();
   printf("*************************************\n");
   return task;
}
