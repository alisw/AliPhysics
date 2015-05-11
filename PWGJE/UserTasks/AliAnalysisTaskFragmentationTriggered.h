#ifndef _ALIANALYSISTASKFRAGMENTATIONTRIGGERED_
#define _ALIANALYSISTASKFRAGMENTATIONTRIGGERED_


class TList;
class TH1;
class TH2;
class TH3;
class AliESDtrackCuts;
class TFormula;
class AliAnalysisUtils;
class TObjArray;
class TRandom3;

class AliAnalysisTaskFragmentationTriggered : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskFragmentationTriggered(const char *name = "fragmentationtriggered"); //default constructor
  
  //AliAnalysisTaskFragmentationTriggered(const AliAnalysisTaskFragmentationTriggered&); //copy constructor
  //AliAnalysisTaskFragmentationTriggered& operator=(const AliAnalysisTaskFragmentationTriggered&); // = operator
  
  virtual ~AliAnalysisTaskFragmentationTriggered(); //destructor

  //  virtual void   Initialize();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(const Option_t *option);

  void SetTrigger(Int_t trg) {
    if (trg == 1) {
      kmbtrigger = kFALSE;
      khjttrigger = kTRUE;
    }
    else {
      kmbtrigger = kTRUE;
      khjttrigger = kFALSE;
    }
  }

  void SetLeading(Double_t value){
    if (value>fLeading) fLeading = value;

  }
  
  double GetLeading(void){
    return fLeading;
  }

  void ClearLeading(void){
    fLeading = 0;
  }
  

protected:
  //  AliAnalysisTaskFragmentationTriggered(const AliAnalysisTaskFragmentationTriggered&);

  AliAnalysisTaskFragmentationTriggered(const AliAnalysisTaskFragmentationTriggered&); //copy constructor 
  AliAnalysisTaskFragmentationTriggered& operator=(const AliAnalysisTaskFragmentationTriggered&); // = operator 

  virtual void Initialize();

  //  TString * trigger;
  Double_t fLeading;
  bool kmbtrigger;
  bool khjttrigger;
  //MC Truth objects
  TObjArray *fMCRC;
  TObjArray *fMCjetstracks;
  TObjArray *fMCjets;
  TObjArray *fMCtracks_pcp;
  TObjArray *fMCtracks_pcm;
  //reco objects
  TObjArray  *fRC;
  TObjArray  *fJets;
  TObjArray  *fJetsTracks;
  TObjArray *fTracks_pcp;
  TObjArray *fTracks_pcm;
  TObjArray *fTracks_pcp_LJCM;
  TObjArray *fTracks_LJCM;
  AliESDtrackCuts *fCuts;
  AliESDtrackCuts *fCuts_noSPD;
  AliESDtrackCuts *fCuts_noITS;
  AliAnalysisUtils *fAnaUtils; 
  TObjArray  *fMCTruthParticles;
  TRandom3   *fRandomGen; //!  dont stream fRandomGen 
 
  TList *fOutputList;
  //MC Histos
  TH1   *fHistMCTruthPtHard;
  TH1   *fHistMCTruthPtHardTrials;
  TH1   *fHistMCTruthNTracks;
  TH1   *fHistMCTruthTrackPt;
  TH1   *fHistMCTruthTrackPhi;
  TH1   *fHistMCTruthTrackEta;
  TH1   *fMCTruthRecTrackPtGen;
  TH1   *fMCTruthRecTrackSecPtGen;
  TH1   *fMCTruthRecTrackSecPtRec;
  TH2   *fHistMCTruthTrackPtRatio;
  TH2   *fHistMCTruthTrackPtRatiowSPD;
  TH2   *fHistMCTruthTrackPtRationoSPD;
  TH1   *fMCTruthnJets;
  TH1   *fMCTruthAntiktHistPt;
  TH2   *fMCTruthAntiktHistAreavsPt;
  TH2   *fMCTruthAntiktHistAreaoverCirclevsPt;
  TH2   *fMCTruthAntiktHistxipt;
  TH2   *fMCTruthAntiktHistxipt_pcp;
  TH2   *fMCTruthAntiktHistxipt_pcm;
  TH2   *fMCTruthAntiktHistxipt_pc;
  TH1   *fMCTruthAntiktHistPtLeading;
  TH1   *fMCTruthAntiktHistEtaLeading;
  TH1   *fMCTruthAntiktHistPhiLeading;
  TH2   *fMCTruthAntiktHistxiptLeading;
  TH2   *fMCTruthAntiktHistxiptLeading_pcp;
  TH1   *fMCTruthRCntry;
  TH1   *fMCTruthAntiktHistPhiRC;
  TH1   *fMCTruthAntiktHistEtaRC;
  TH1   *fMCTruthRCnTracksOverlap;
  TH1   *fMCTruthHistRCPt;
  TH1   *fMCTruthHistRCnTracks;
  TH1   *fMCTruthHistRCWithOverlapPt;
  TH1   *fMCTruthHistRCWithOverlapnTracks;
  TH1   *fMCTruthHistRCTrackEta;
  TH1   *fMCTruthHistRCTrackPhi;
  TH1   *fMCTruthHistRCTrackEtaDiscarded;
  TH1   *fMCTruthHistRCTrackPhiDiscarded;
  TH2   *fMCTruthHistRCxipt;
  TH2   *fMCTruthHistRCxiptLeading;
  TH2   *fMCTruthHistRCxiptLeadingDiscarded;
  //rho histos
  TH1   *fMCTruthHistRho;
  TH2   *fMCTruthHistRhovsLJPt;
  TH1   *fMCTruthHistRhoSubPt;
  TH2   *fMCTruthHistRhoSubAreaOverCirclevsPt;
  TH2   *fMCTruthHistRhoSubPhiTracksJetvsPt;
  TH2   *fMCTruthHistRhoSubEtaTracksJetvsPt;
  TH2   *fMCTruthHistRhoSubPhiTracksPCvsPt;
  TH2   *fMCTruthHistRhoSubEtaTracksPCvsPt;
  TH2   *fMCTruthHistRhoSubXiPt;
  TH2   *fMCTruthHistRhoSubXiPCPt;



  //track histos
  TH1   *fHistSumw2Debug;
  TH1   *fHistEventDebug;
  TH1   *fHistCuts;
  TH1   *fHistTrackCuts;
  TH1   *fHistTrackCutsDep;
  TH1   *fHistPt;
  TH1   *fHistnTracks;
  TH1   *fHistEta;
  TH1   *fHistPhi;
  TH2   *fHistEtaPhi;
  TH1   *fHistPhiHybridUnconstrained;
  TH1   *fHistPhiHybridConstrainedwITS;
  TH1   *fHistPhiHybridConstrainedwoITS;
  TH1   *fHistITSNcls;
  TH1   *fHistITShitsinlayer;
  TH2   *fHistSPDNclsvstracklets;
  TH1   *fHistITSchi2overhits;
  TH2   *fHistMulthybridvsglobal; 
  TH1   *fHistTPCNcls;
  TH1   *fHistTPCchi2;
  TH1   *fHistTPCchi2overtpchits; 
  //  TH3   *fHistVertex;
  TH1   *fHistImpactParameterglobalxy;
  TH1   *fHistImpactParameterglobalz;
  TH1   *fHistImpactParameternoSPDxy;
  TH1   *fHistImpactParameternoSPDz;
  TH1   *fHistImpactParameternoITSxy;
  TH1   *fHistImpactParameternoITSz;
  /*
  TH1   *fHistDCAxyglobal_esd;
  TH1   *fHistDCAzglobal_esd;
  TH1   *fHistDCAxywoSPD_esd;
  TH1   *fHistDCAzwoSPD_esd;
  TH1   *fHistDCAxywoITS_esd;
  TH1   *fHistDCAzwoITS_esd;
  */
  TH2   *fHistVertexXY;
  TH1   *fHistVertexZ;
  TH1   *fHistContribtoVtx;
  /*
  TH1   *fHistDCAxyglobal;
  TH1   *fHistDCAzglobal;
  TH2   *fHistDCAzvsVTXglobal;
  TH1   *fHistDCAxywoSPD;
  TH1   *fHistDCAzwoSPD; 
  TH2   *fHistDCAzvsVTXwoSPD;
  TH1   *fHistDCAxywoITS;
  TH1   *fHistDCAzwoITS;
  TH1   *fHistDCAzvsVTXwoITS;
  TH2   *fHistXYAtDCAwoSPD;
  TH1   *fHistZAtDCAwoSPD;
  TH2   *fHistXYAtDCAwoITS;
  TH1   *fHistZAtDCAwoITS;  
  */
  //jet histos
  TH1   *fHistnJets;
  TH1   *fAntiktHistPt;
  TH1   *fAntiktHistPtwTRD;
  TH1   *fAntiktHistPtwoTRD;
  TH1   *fAntiktHistEta;
  TH2   *fAntiktHistEtavsPt;
  TH1   *fAntiktHistPhi;
  TH2   *fAntiktHistEtaPhi;
  TH2   *fAntiktHistAreavsPt;
  TH2   *fAntiktHistAreaoverCirclevsPt;
  TH1   *fAntiktHistnTracks;
  TH1   *fAntiktHistz;
  TH2   *fAntiktHistzpt;
  TH1   *fAntiktHistxi;
  TH2   *fAntiktHistxipt;
  TH2   *fAntiktHistxipt_pcp;
  TH2   *fAntiktHistxipt_pcm;
  TH2   *fAntiktHistxipt_pc;
  TH2   *fAntiktHistxipt_RC;
  TH2   *fAntiktHistxiptwTRD;
  TH2   *fAntiktHistxiptwoTRD;

  //Random Cone histos
  TH1   *fHistRCntry;
  TH1   *fHistPhiRC;
  TH1   *fHistEtaRC;
  TH1   *fHistRCnTracksOverlap;
  TH1   *fHistRCPt;
  TH1   *fHistRCnTracks;
  TH2   *fHistRCnTracksvsLJPt;
  TH2   *fHistRCPtvsLJPt;

  //leading jet histos
  //  TH1   *fHistnJetsLeading;
  TH1   *fAntiktHistPtLeading;
  TH2   *fAntiktHistPtLeadingvsHybridMult;
  TH1   *fAntiktHistEtaLeading;
  TH1   *fAntiktHistPhiLeading;
  TH2   *fAntiktHistEtaPhiLeading;
  TH1   *fAntiktHistAreaLeading;
  TH1   *fAntiktHistAreaoverCircleLeading;
  TH2   *fAntiktHistAreaLeadingvsPt;
  TH2   *fAntiktHistAreaoverCircleLeadingvsPt;
  TH2   *fAntiktHistAreaoverCircleLeadingvsNtracks;
  TH1   *fAntiktHistnTracksLeading;
  TH2   *fAntiktHistnTracksLeadingvsPt;
  TH1   *fAntiktHistTrackPtLeading;
  TH1   *fAntiktHistTrackPhiLeading;
  TH1   *fAntiktHistTrackEtaLeading;
  TH1   *fAntiktHistzLeading;
  TH2   *fAntiktHistzptLeading;
  TH1   *fAntiktHistxiLeading;
  TH2   *fAntiktHistxiptLeading;
  TH2   *fAntiktHistxitimesptptLeading;
 
  //Perpendicular Cone Histos:
  //Plus histos:
  TH1   *fLJPCPTrackPt;
  TH1   *fLJPCPTrackEta;
  TH1   *fLJPCPTrackPhi;
  TH2   *fLJPCPTrackEtaPhi;
  TH1   *fLJPCPnTracks;
  TH2   *fLJPCPnTracksPhi;
  TH1   *fLJPCPPt;
  TH2   *fLJPCPPtvsLJPt;
  TH1   *fLJPCPEta;
  TH1   *fLJPCPPhi;
  TH2   *fLJPCPzpt;
  TH2   *fLJPCPxipt;



  //Minus histos:
  TH1   *fLJPCMTrackPt;
  TH1   *fLJPCMTrackEta;
  TH1   *fLJPCMTrackPhi;
  TH2   *fLJPCMTrackEtaPhi;
  TH1   *fLJPCMnTracks;
  TH2   *fLJPCMnTracksPhi;
  TH1   *fLJPCMPt;
  TH2   *fLJPCMPtvsLJPt;
  TH1   *fLJPCMEta;
  TH1   *fLJPCMPhi;
  TH2   *fLJPCMzpt;
  TH2   *fLJPCMxipt;

  //both cones
  TH1   *fLJPCPMMPt;
  TH1   *fLJPCPt;
  TH2   *fLJPCPtvsLJPt;

  //perp cone subtracted
  TH1   *fAntiktHistPtLeading_sub_pcp;
  TH2   *fAntiktHistxiptLeading_sub_pcp;
  TH1   *fAntiktHistTrackPtLeading_sub_pcp;

  TH2   *fAntiktHistxiptLeading_pcp_sub_pcp;
  TH1   *fAntiktHistTrackPtLeading_pcp_sub_pcp;

  //leading jet cone method histos:
  TH2   *fAntiktHistxiptLeading_pcp_LJCM;
  TH1   *fAntiktHistPtLeading_LJCM;
  TH2   *fAntiktHistxiptLeading_LJCM;
  TH1   *fAntiktHistTrackPtLeading_LJCM;

  //leading track cut
  TH2   *fAntiktHistxiptLeading_LTC;
  TH2   *fAntiktHistxiptLeading_LTC_pcp;
  TH1   *fAntiktHistPtLeading_LTC;
 
 // Event by event UE pT subtraction

  TH1   *fHistRho;
  TH2   *fHistRhovsLJPt;
  TH1   *fHistRhoSubPt;
  TH2   *fHistRhoSubAreaOverCirclevsPt;
  TH2   *fHistRhoSubPhiTracksJetvsPt;
  TH2   *fHistRhoSubEtaTracksJetvsPt;
  TH2   *fHistRhoSubPhiTracksPCvsPt;
  TH2   *fHistRhoSubEtaTracksPCvsPt;
  TH2   *fHistRhoSubXiPt;
  TH2   *fHistRhoSubXiPCPt;


  // Event by event RC UE pT subtraction                                                                                                                                               

  TH1   *fHistRhoRC;
  TH2   *fHistRhoRCvsLJPt;
  TH1   *fHistRhoRCSubPt;
  TH2   *fHistRhoRCSubAreaOverCirclevsPt;
  TH2   *fHistRhoRCSubPhiTracksJetvsPt;
  TH2   *fHistRhoRCSubEtaTracksJetvsPt;
  TH2   *fHistRhoRCSubPhiTracksPCvsPt;
  TH2   *fHistRhoRCSubEtaTracksPCvsPt;
  TH2   *fHistRhoRCSubXiPt;
  TH2   *fHistRhoRCSubXiPCPt;


  //delta pt Histos
  TH1   *fHistDeltaPt;
  TH2   *fHistDeltaPtvsLJPt;


  ClassDef(AliAnalysisTaskFragmentationTriggered, 1);
};

#endif
