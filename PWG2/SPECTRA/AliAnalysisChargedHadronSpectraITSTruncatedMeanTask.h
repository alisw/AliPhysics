//Class to extract data to do ITS+TPC global Spectra
//Autor Marek Chojnacki
//emali Marek.Chojnacki@cern.ch



#ifndef ALIANALYSISCHARGEDHADRONSPECTRAITSTRUNCATEDMEANTASK_H
#define  ALIANALYSISCHARGEDHADRONSPECTRAITSTRUNCATEDMEANTASK_H
//#include <fstream>
class TH1F;
class TH2F;
class TH3F;
class AliESDEvent;
class AliESDtrack;
class  AliESDtrackCuts;
class AliESDpidCuts;
class AliESDpid;
class TGraph;
class AliStack;
#include "AliAnalysisTaskSE.h"
//#include "AliESDtrackCuts.h"
#include "AliESDpidCuts.h"
#include "AliESDpid.h"




class AliAnalysisChargedHadronSpectraITSTruncatedMeanTask : public AliAnalysisTaskSE {
 public:
  AliAnalysisChargedHadronSpectraITSTruncatedMeanTask(const char *name = "AliAnalysisChargedHadronSpectraITSTruncatedMeanTask");
  virtual ~AliAnalysisChargedHadronSpectraITSTruncatedMeanTask() {}
  
  //virtual void   ConnectInputData(Option_t *);
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *); 
  virtual void   LocalInit();
  AliESDtrackCuts*  GetAliESDtrackCuts() const {return fCuts;}
  void SetMCOn(){fMC=kTRUE;
   fESDpid->GetTPCResponse().SetBetheBlochParameters(2.15898e+00/50.,1.75295e+01,3.40030e-09,1.96178e+00,3.91720e+00);}
  void SetAliESDtrackCuts(AliESDtrackCuts* const cuts ){fCuts=cuts;/*flist->Add(fCuts);*/}
  void SetFunctionParam(Double_t* const par);
   void SetMultiplicityCut(Int_t low, Int_t up);
 void SetSPDMethodCut(){fSPD=kTRUE;}
  void SetCorrectSDD(){fCorrectSDD=kTRUE;}
   void SetCorrectSSD(){fCorrectSSD=kTRUE;}
  void SetYcut(Float_t value){fYCut=TMath::Abs(value);}
  void Setsigmacut(Float_t value){fsigmacut=TMath::Abs(value);}
  void SetNsigmaDCAcut(Float_t sigmaxy,Float_t sigmaz){fnsigmaxy=sigmaxy;fnsigmaz=sigmaz;}    
  void SetChargeCut(Float_t chargeCut){fchargeCut=TMath::Abs(chargeCut)>50.0?50.0:TMath::Abs(chargeCut);}
     void SetTPCPIDCUT(AliESDpidCuts* const cuts){fTPCPIDCUT=cuts;}
     void SetWeights(TGraph* const setK0weight, TGraph* const  setlambdaweight,TGraph* const setAntilambdaweight){fK0weight=setK0weight;flambdaweight=setlambdaweight;fAntilambdaweight=setAntilambdaweight;}
     void SetDCA2010();
     void SetHImode(){fHIsettings=kTRUE;}
     void SetCentralityCut(Float_t low, Float_t up); 
     void SetDoVertexrescuts(){fdovertexrescuts=kTRUE;}
 private:
 
 
 AliESDEvent *fESD;    //ESD object    
  AliESDtrackCuts *fCuts;//cuts 
  AliESDtrackCuts *fCutsMul;//cuts for multiplicty 
  
  
  
  Bool_t fMC;//if TRUE use MC 
  Int_t fLowMultiplicity;//low Multiplicity cut
  Int_t fUpMultiplicity;//up Multiplicity cut
  Float_t fLowCentrality;//low Centrality cut
  Float_t fUpCentrality;//up  Centrality cut
  Bool_t fSPD;//use spd2 as mulestimator 
  
  Float_t fYCut;//cut in y
  Float_t fsigmacut;//cut in sigma in n-sigma method
  Float_t fnsigmaxy; //cut in sigma on xy dca
  Float_t fnsigmaz;//cut in sigma on Z dca
  Float_t fdcaxypar[3];//parameters for DCAxy cut 
  Float_t fdcazpar[4];//parameters for DCAz cut 
  
  Float_t fchargeCut;//cut for the low charges
  
  
 Bool_t fCorrectSDD;//In LHC10a3 in some runs dE in SDD had to scaled to SSD flag if this should be done 
  Bool_t fCorrectSSD;//this same but for dE SSS
  
  Bool_t fHIsettings;//speciall settings fot HI mode
  Bool_t fdovertexrescuts;// check on Vmc-VESD 
  
  
  
  TGraph* fK0weight ;//weight for pions comming from K0shorts
  TGraph* flambdaweight ;//weight for protons comming from lambdas
   TGraph* fAntilambdaweight ;//weight for antiprotons comming from antilambdas
  
   
   
   
TH1F *fHistStats; //histogram with statistic of events
TH1F* fHistZVertexBeforeCut; //Z of vertex before cut 
TH1F* fHistZVertexAfterCut; //Z of vertex after cut
TH2F* fHistXYVertexBeforeCut; //XY of vertex before cut 
TH2F* fHistXYVertexAfterCut; //XY of vertex after cut


TH2F* fHistPhiPtBeforeCuts;//phi pt before cuts 
TH2F* fHistPhiPtAfterCuts;//phi pt after cuts 
TH2F* fHistEtaPtBeforeCuts;//eta pt before cuts 
TH2F* fHistEtaPtAfterCuts;//eta pt after cuts 

TH2F* fHistDCABeforeCuts;//dca hist before cuts 
TH2F* fHistDCAAfterCuts;//dca hist after cuts 

TH2F* fHistPminusTPCinPAfterCuts;//differnece between global momentum at primary vetrex and tpc standalone momentum at primry vertex 
TH2F* fHistPminusTPCinPglobalAfterCuts;//differnece between global momentum at primary vetrex and global momentum at  the inner  wall of the TPC taken from global tracking

//positive
TH2F* fHistMydEPpositive;//dE in its as function of global p at p.v.
TH2F* fHistMydETPCinPpositive;//dE in its as function of TPC p at p.v.
TH2F* fHistMydETPCinPglobalpositive;//dE in its as function of p  at  the inner  wall of the TPC taken from global tracking
//negative
TH2F* fHistMydEPnegative;//dE in its as function of global p at p.v.
TH2F* fHistMydETPCinPnegative;//dE in its as function of TPC p at p.v.
TH2F* fHistMydETPCinPglobalnegative;//dE in its as function of p  at  the inner  wall of the TPC taken from global tracking

//dE   as function of global p at p.v.
TH2F* fHistL3dEP;// SDD1 
TH2F* fHistL4dEP;// SDD2
TH2F* fHistL5dEP;//SSD1
TH2F* fHistL6dEP; //SSD2

//dE in its as function of TPC p at p.v.
TH2F* fHistL3dETPCinP;// SDD1 
TH2F* fHistL4dETPCinP;// SDD2
TH2F* fHistL5dETPCinP;//SSD1
TH2F* fHistL6dETPCinP;//SSD2




TH2F* fHistwhichhasmin;// ITS  layer with minimal charged
TH1F* fHistMysignalminusESD;// My signal minus ESD


//log dE-logdEfit as function of  global p at p.v. for
TH2F* fHistminsignalifPionP;//pions
TH2F* fHistminsignalifKaonP;//kaons
TH2F* fHistminsignalifProtonP;//protons

TH2F* fHistminsignalifAntiPionP;//antipions
TH2F* fHistminsignalifAntiKaonP;//antikaons
TH2F* fHistminsignalifAntiProtonP;//antiprotons


//DCA histograms for clean particles with after dca cut
TH3F* fDCAXYZforcleanPions;//pions
TH3F* fDCAXYZforcleanAntiPions;//antipion
TH3F* fDCAXYZforcleanProtons;//kaons
TH3F* fDCAXYZforcleanAntiProtons;//antikaons

//DCA histograms for clean particles with before dca cut
TH3F* fDCAXYZOpenforcleanPions;//pions
TH3F* fDCAXYZOpenforcleanAntiPions;//antipions
TH3F* fDCAXYZOpenforcleanProtons;//kaons
TH3F* fDCAXYZOpenforcleanAntiProtons;//antikaons

//pt distibution of track fullfilling some cuts 0-pions, 1-kaons,2-protons
TH2F* fHistNtrackwithstandardcuts;//TPC cuts
TH2F* fHistNtrackwithITSPIDcuts;//TPC cuts + ITS pid cuts

TH2F* fHistSignalinTPCKaonforstandardcuts;//TPC signal for Kaons tpc cuts tracks
TH2F* fHistSignalinTPCKaonforITSPIDcuts;//TPC signal for Kaons tpc+itspid cuts tracks

TH2F* fHistSignalinTPCAntiKaonforstandardcuts;//TPC signal for AntiKaons tpc cuts tracks
TH2F* fHistSignalinTPCAntiKaonforITSPIDcuts; //TPC signal for AntiKaons tpc+itspid cuts tracks


TH2F* fHistSignalinTPCProtonforstandardcuts; //TPC signal for Protons tpc cuts tracks
TH2F* fHistSignalinTPCProtonforITSPIDcuts;//TPC signal for Protons tpc+itspid cuts tracks

TH2F* fHistSignalinTPCAntiProtonforstandardcuts;//TPC signal for AntiProtons tpc cuts tracks
TH2F* fHistSignalinTPCAntiProtonforITSPIDcuts; //TPC signal for AntiProtons tpc+itspid cuts tracks


//Multiplicity histos
TH1F* fHistStandartMul;//number from AliESDtrackCuts::GetReferenceMultiplicity
TH1F* fHistMytrackMul;//number of my tracks
TH2F* fHistStandartMulvSPD2;//number from AliESDtrackCuts::GetReferenceMultiplicity v SPD2

//log dE-logdEfit as function of  global p at p.v. for primary tracks 
TH2F* fHistminsignalifPionPPrimary; //pions
TH2F* fHistminsignalifKaonPPrimary;//kaons
TH2F* fHistminsignalifProtonPPrimary;//protons
TH2F* fHistminsignalifProtonPPrimaryfake;//fake protons

TH2F* fHistminsignalifAntiPionPPrimary;//antipions
TH2F* fHistminsignalifAntiKaonPPrimary;//antikaons
TH2F* fHistminsignalifAntiProtonPPrimary;//antiprotons
TH2F* fHistminsignalifAntiProtonPPrimaryfake;//antiprotonsfake


//log dE-logdEfit as function of  global p at p.v. for secondary and other tracks tracks 
TH2F* fHistminsignalifPionPSecondary;//pions
TH2F* fHistminsignalifKaonPSecondary;//kaon
TH2F* fHistminsignalifProtonPSecondaryWD;//protons comming from weak decays 
TH2F* fHistminsignalifProtonPSecondaryHI;//protons comming from material 
TH2F* fHistminsignalifProtonPSecondaryRest;//rest contamination 

TH2F* fHistminsignalifProtonPSecondaryWDfake;//protons fakes comming from weak decays 
TH2F* fHistminsignalifProtonPSecondaryHIfake;//protons fakes comming from material 

TH2F* fHistminsignalifAntiPionPSecondary;//antipions
TH2F* fHistminsignalifAntiKaonPSecondary;//antikaon

TH2F* fHistminsignalifAntiProtonPSecondaryWD;//antiprotons comming from weak decays
TH2F* fHistminsignalifAntiProtonPSecondaryHI;//antiprotons fakes comming from material 
TH2F* fHistminsignalifAntiProtonPSecondaryRest;//rest contamination 

TH2F* fHistminsignalifAntiProtonPSecondaryWDfake;//antiprotons fakes comming from weak decays 
TH2F* fHistminsignalifAntiProtonPSecondaryHIfake;//antiprotons fakes comming from material 
 

TH2F* fHistminsignalifMuEPositiveP;//mu+  positronium for pions  
TH2F* fHistminsignalifMuENegativeP;//mu-  electrons for antipions 

 TH2F* fHistminsignalifPionPrimaryfake;//fake pions primary
  TH2F* fHistminsignalifKaonPrimaryfake;//fake kaons primary
 
 TH2F* fHistminsignalifAntiPionPrimaryfake;//fake antipions primary
 TH2F* fHistminsignalifAntiKaonPrimaryfake;//fake antikaons primary
 
 
  TH2F* fHistminsignalifPionSecondaryfake;//fake pions
  TH2F* fHistminsignalifKaonSecondaryfake;//fake kaons
 
 TH2F* fHistminsignalifAntiPionSecondaryfake;//fake antipions
 TH2F* fHistminsignalifAntiKaonSecondaryfake;//fake antikaons
 

 //MC particles from Events passing ESD event cuts only pt
TH1F* fHistminsignalifPionPMCPrimary;//Pions 
TH1F* fHistminsignalifKaonPMCPrimary;//Kaons
TH1F* fHistminsignalifProtonPMCPrimary;//Protons

TH1F* fHistminsignalifAntiPionPMCPrimary;//AntiPions
TH1F* fHistminsignalifAntiKaonPMCPrimary;//AntiKaons
TH1F* fHistminsignalifAntiProtonPMCPrimary;//AntiProtons


//MC particles from all MC events 
TH1F* fHistminsignalifPionPMCPrimaryBeforeEventCuts;//Pions
TH1F* fHistminsignalifKaonPMCPrimaryBeforeEventCuts;//Kaons
TH1F* fHistminsignalifProtonPMCPrimaryBeforeEventCuts;//Protons

TH1F* fHistminsignalifAntiPionPMCPrimaryBeforeEventCuts;//AntiPions
TH1F* fHistminsignalifAntiKaonPMCPrimaryBeforeEventCuts;//AntiKaons
TH1F* fHistminsignalifAntiProtonPMCPrimaryBeforeEventCuts;//AntiProtons

//MC particles from all MC events  with good vertex in z 
TH1F* fHistminsignalifPionPMCPrimaryBeforeEventCutswithgoodZvertex;//Pions
TH1F* fHistminsignalifKaonPMCPrimaryBeforeEventCutswithgoodZvertex;//Kaons
TH1F* fHistminsignalifProtonPMCPrimaryBeforeEventCutswithgoodZvertex;//Protons

TH1F* fHistminsignalifAntiPionPMCPrimaryBeforeEventCutswithgoodZvertex;//AntiPions
TH1F* fHistminsignalifAntiKaonPMCPrimaryBeforeEventCutswithgoodZvertex;//AntiKaons
TH1F* fHistminsignalifAntiProtonPMCPrimaryBeforeEventCutswithgoodZvertex;//AntiProtons


//MC particles from  MC events  which ESD event go trought physics selection and has vertex  
TH1F* fHistminsignalifPionPMCPrimaryAfterEventCutsBeforeVertexZ;//Pions
TH1F* fHistminsignalifKaonPMCPrimaryAfterEventCutsBeforeVertexZ;//Kaons
TH1F* fHistminsignalifProtonPMCPrimaryAfterEventCutsBeforeVertexZ;//Protons

TH1F* fHistminsignalifAntiPionPMCPrimaryAfterEventCutsBeforeVertexZ;//AntiPions
TH1F* fHistminsignalifAntiKaonPMCPrimaryAfterEventCutsBeforeVertexZ;//AntiKaons
TH1F* fHistminsignalifAntiProtonPMCPrimaryAfterEventCutsBeforeVertexZ;//AntiProtons


//DCA distributions for different parctiles after dca cuts 

TH3F* fDCAXYZforcleanPionsMCPrimary;//primary pions 
TH3F* fDCAXYZforcleanAntiPionsMCPrimary;//primary antipions 
TH3F* fDCAXYZforcleanProtonsMCPrimary;//primary protons
TH3F* fDCAXYZforcleanAntiProtonsMCPrimary;//primary antiprotons

//Secondrary Pions weak deacy
TH3F* fDCAXYZforcleanPionsWD;//pions 
TH3F* fDCAXYZforcleanAntiPionsWD;// antipions 

//Secondrary Protons weak deacy + fakes
TH3F* fDCAXYZforcleanProtonsWD;//protons
TH3F* fDCAXYZforcleanAntiProtonsWD;//antiprotons

//Secondrary Pions Hadronic
TH3F* fDCAXYZforcleanPionsHI;//pions
TH3F* fDCAXYZforcleanAntiPionsHI;//antipions

//Secondrary Protons Hadronic+fakes
TH3F* fDCAXYZforcleanProtonsHI;//pions
TH3F* fDCAXYZforcleanAntiProtonsHI;//antipions

//Secondrary Pions mu el
TH3F* fDCAXYZforcleanPionsMEPrimary;//posvitive
TH3F* fDCAXYZforcleanAntiPionsMEPrimary;//negative
TH3F* fDCAXYZforcleanPionsMESecondary;//posvitive
TH3F* fDCAXYZforcleanAntiPionsMESecondary;//negative

//Secondrary Pions rest source
TH3F* fDCAXYZforcleanPionsR;//positive
TH3F* fDCAXYZforcleanAntiPionsR;//negative

//Secondrary Protons rest
TH3F* fDCAXYZforcleanProtonsR;//positive
TH3F* fDCAXYZforcleanAntiProtonsR;//negative


//DCA distributions for different parctiles before dca cuts 

TH3F* fDCAXYZOpenforcleanPionsMCPrimary;//primary pions 
TH3F* fDCAXYZOpenforcleanAntiPionsMCPrimary;//primary antipions 
TH3F* fDCAXYZOpenforcleanProtonsMCPrimary;//primary protons
TH3F* fDCAXYZOpenforcleanAntiProtonsMCPrimary;//primary antiprotons

//Secondrary Pions weak deacy
TH3F* fDCAXYZOpenforcleanPionsWD;//pions 
TH3F* fDCAXYZOpenforcleanAntiPionsWD;// antipions 

//Secondrary Protons weak deacy + fakes
TH3F* fDCAXYZOpenforcleanProtonsWD;//protons
TH3F* fDCAXYZOpenforcleanAntiProtonsWD;//antiprotons

//Secondrary Pions Hadronic
TH3F* fDCAXYZOpenforcleanPionsHI;//pions
TH3F* fDCAXYZOpenforcleanAntiPionsHI;//antipions

//Secondrary Protons Hadronic+fakes
TH3F* fDCAXYZOpenforcleanProtonsHI;//pions
TH3F* fDCAXYZOpenforcleanAntiProtonsHI;//antipions

//Secondrary Pions mu el
TH3F* fDCAXYZOpenforcleanPionsMEPrimary;//posvitive
TH3F* fDCAXYZOpenforcleanAntiPionsMEPrimary;//negative
TH3F* fDCAXYZOpenforcleanPionsMESecondary;//posvitive
TH3F* fDCAXYZOpenforcleanAntiPionsMESecondary;//negative

//Secondrary Pions rest source
TH3F* fDCAXYZOpenforcleanPionsR;//positive
TH3F* fDCAXYZOpenforcleanAntiPionsR;//negative

//Secondrary Protons rest
TH3F* fDCAXYZOpenforcleanProtonsR;//positive
TH3F* fDCAXYZOpenforcleanAntiProtonsR;//negative




//Electron Muon  source procces
TH2F* fElectronsource; //e+
TH2F* fAntiElectronsource;//e-

TH2F* fMuonsource;//mu+
TH2F* fAntiMuonsource;//mu-

//N tpc clusters for 
TH2F* fPionNTPCClusters; //pions tracks
TH2F* fAntiPionNTPCClusters;//antipions tracks 

TH2F* fKaonNTPCClusters; //Kaons tracks
TH2F* fAntiKaonNTPCClusters;//antiKaons tracks 

TH2F* fProtonNTPCClusters; //Protons tracks
TH2F* fAntiProtonNTPCClusters;//antiProtons tracks

TH2F* fPionchi2; //pions tracks
TH2F* fAntiPionchi2;//antipions tracks 

TH2F* fKaonchi2; //Kaons tracks
TH2F* fAntiKaonchi2;//antiKaons tracks 

TH2F* fProtonchi2; //Protons tracks
TH2F* fAntiProtonchi2;//antiProtons tracks
 

TH2F* fTracksCutmonitoring;// Number of tracks as fun of pt on each step of selection 
TH3F* fParticlesCutmonitoring;//Number of as particles as fun of pt on each step of selection x 0-pion 1-kaon,2-proton,3-antipion,4-antikaon,5-antiproton
TH3F* fVertexshift; //shift of the vertex due to reconstruction

TH3F* fPtESDminusPtMCvPtESDafterallcuts;//ptESD -ptMC v ptESD after all cuts 
TH3F* fPtESDminusPtMCvPtESDafterTPCcuts;//ptESD - ptMC v ptESD after TPC cuts (refit,chi2,nclus);

TH3F* fMulESDMulMCVz;//Multiplicty ESD Multiplicty MC Vrt Z 



//TPC pid objects 
AliESDpidCuts* fTPCPIDCUT;//cut
AliESDpid* fESDpid; // global thing

TH1F* fPrimaryElectronsMother; //name says all 





 TList *flist;//output list
	
 AliAnalysisChargedHadronSpectraITSTruncatedMeanTask(const AliAnalysisChargedHadronSpectraITSTruncatedMeanTask&); // not implemented
AliAnalysisChargedHadronSpectraITSTruncatedMeanTask& operator=(const AliAnalysisChargedHadronSpectraITSTruncatedMeanTask&); // not implemented


 Float_t MyITSsignalusing4points(Double_t* const) const;
  Float_t MyITSsignalusing3points(Double_t* const) const;
void CorrectSDD(Double_t *tmpQESD) const;
void CorrectSSD(Double_t *tmpQESD) const;
    Bool_t SelectOnImpPar(AliESDtrack* const t) const;
    Float_t GetWeight(Int_t type,AliStack* const stack) const;
    
    
    
 ClassDef(AliAnalysisChargedHadronSpectraITSTruncatedMeanTask, 2); // example of analysis
};

#endif
