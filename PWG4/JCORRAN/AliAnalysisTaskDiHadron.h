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

//2- and 3-particle trigger particle correlation analysis
//Author: Jason Glyndwr Ulery, ulery@uni-frankfurt.de
//version 3.4, last revised 2010/08/15

#ifndef ALIANALYSISTASKDIHADRON_H
#define ALIANALYSISTASKDIHADRON_H


class TF1;
class TH1F;
class TH2F;
class TH3F;
class AliESDEvent;

#include "AliAnalysisTask.h"


class AliAnalysisTaskDiHadron : public AliAnalysisTask{
  //Default constructor
 public:


  AliAnalysisTaskDiHadron(const char *name="AliAnalysisTaskDiHadron");

  virtual ~AliAnalysisTaskDiHadron() {}
  //virtual ~AliAnalysisTaskDiHadron();

 

  //AliAnalysisTask Functions
 
  virtual void ConnectInputData(Option_t *option);
  virtual void CreateOutputObjects();
  virtual void Exec(Option_t *option);
  virtual void Terminate(Option_t *);
 

  void SetEfficiencies(Float_t EffFitPt, const TF1 *FitLow, const TF1 *FitHigh, Int_t NFitLowParam, Int_t NFitHighParam, Float_t *FitLowParam, Float_t *FitHighParam);
  void SetBins(Int_t nBinPhi, Int_t nBinEta, Int_t nBinPhiEtaPhi, Int_t nBinPhiEtaEta, Int_t nBinPhi3, Int_t nBinEta3, Float_t dPhiMin, Float_t dPhiMax, Int_t NTPtBins, Int_t NMixBins, Int_t NCentBins,Int_t NAPtBins, Int_t NAPt3Bins, Int_t NVertexBins, Int_t NXEBin,Float_t *PtTrigArray, Float_t *PtAssocArray, Float_t *PtAssoc3Array1, Float_t *PtAssoc3Array2, Int_t *CentArrayMin, Int_t *CentArrayMax, Float_t *XEArray);
  void SetOptions(Int_t fEfficiencyCorr, Int_t fDEBUG,Int_t fMCHistos);
  void SetCuts(Int_t MinClutersTPC, Float_t MinClusterRatio, Float_t MaxTPCchi2, Int_t MinClustersITS, Float_t EtaCut, Float_t TrigEtaCut, Float_t NearPhiCut, Float_t XECut, Float_t MaxDCA, Float_t MaxDCAXY, Float_t MaxDCAZ, Int_t DCA2D, Int_t TPCRefit, Int_t ITSRefit, Int_t SPDCut, Float_t MinPtAssoc, Float_t MaxPtAssoc, Float_t VzCut, Int_t NIDs, const char *TrigIDArray);
  Int_t CheckVertex(const AliESDEvent *rESD);
  Int_t CheckTrigger(const AliESDEvent *rESD);
  Int_t TrackCuts(const AliESDEvent *rESD, Float_t *rPt, Float_t *rEta, Float_t *rPhi, Short_t *rCharge, Float_t *rEff, Int_t **rPtAssoc3, Int_t *rNPtAssoc3, Int_t *rGoodTracks);
  Int_t TrackCutsMC(AliMCEvent *rMC, Float_t *rPt, Float_t *rEta, Float_t *rPhi, Short_t *rCharge, Float_t *rEff, Int_t **rPtAssoc3, Int_t *rNPtAssoc3, Int_t *rGoodTracks);

  
 private:

//Maximum numbers for creating the arrays
  enum{kNumberOfPtBins=20, 
       kNumberOfCentBins=10, 
       kNumberOfAPtBins=50,
       kNumberOfApt3Bins=50,
       kNumberOfTriggerIDs=10,
       kNumberOfVertexBins=20,
       kNumberOfXeBins=20,
       kNumberOfEventsToMix=100};
  // virtual void SetParameters();
  AliESDEvent *fESD; //ESD data file
  AliMCEvent *fMC;//Monty Carlo event file
  TList *fOutput;//List used to store all the histograms for output
  
   //Cuts
  
  Int_t fMinClustersTPC;//Minimum Clusters in the TPC
  Float_t fMinClusterRatio;//Ratio of MinClusters to Number of Possible Clusters in TPC
  Float_t fMaxTPCchi2;//Max chi^2/ndf for TPC fit
  Int_t fMinClustersITS;//Minimum Clusters in ITS
  Float_t fEtaCut;//Associated particle eta range |eta|<fEtaCut
  Float_t fTrigEtaCut;//Trigger particle eta range |eta|<fTrigEtaCut (should be <=fEtaCut)
  Float_t fNearPhiCut;//Point in DeltaPhi seperating near and away sides
  Float_t fXECut;//Point in DeltaPhi seperating near and away for XE (not finished)
  Float_t fMaxDCA;//Max total DCA (if fDCA2D==0) (in cm)
  Float_t fMaxDCAXY;//Max radial DCA (if fDCA2D==1)
  Float_t fMaxDCAZ;//Max longitudional DCA (if fDCA2D==1)
  Int_t fDCA2D;//1 cut on total DCA, 2 cut seperately on radial and logitudional, 3 pt dependent cut
  Int_t fTPCRefit;//1 require TPCRefit 0 not
  Int_t fITSRefit;//1 require ITSRefit 0 not
  Int_t fSPDCut;//1 require a point in the SPD 0 not
  Float_t fMinPtAssoc;//Minimum pT to look at 
  Float_t fMaxPtAssoc;//Maximum pT to look at  
  Float_t fVzCut;//z-vertex cut |z-vertex|<fVzCut in cm
    Int_t fEfficiencyCorr;//Toggle correcting of efficiencies when filling histograms
    Int_t fDEBUG;//if 1 extra printfs for debugging

    //Binning
    Int_t fnBinPhi;//Number of bins for DeltaPhi histograms
    Int_t fnBinEta;//Number of bins for DeltaEta histograms
    Int_t fnBinPhiEtaPhi;//Number of bins for DeltaPhi-DeltaEta histograms in Phi
    Int_t fnBinPhiEtaEta;//Number of bins for DeltaPhi-DeltaEta histograms in Eta
    Int_t fnBinPhi3;//Number of bins in 3-particle DeltaPhi-DeltaPhi (along each axis)
    Int_t fnBinEta3;//Number of bins in 3-particle DeltaEta-DeltaEta (along each axis)
    Float_t fPi;//Pi=3.14....
    Float_t fdPhiMin;//lowest edge of the first bin in DeltaPhi plots
    Float_t fdPhiMax;//maximum in DeltaPhi plots (should be 2*Pi+fdPhiMin)

    //Parameters for settings
    Int_t fNTPtBins;//Number of trigger pt bins
    Int_t fNMix;//Number of events to mix
    Int_t fNCentBins;//Number of centrality bins (they are allowed to overlap)
    Int_t fNAPtBins;//Number of associated particle bins for 2-particle correlations (overlap not allowed)
    Int_t fNAPt3Bins;//Number of associated particle bins for 3-particle correlations (overlap allowed)
    Int_t fNVertexBins;//Number of bins in z-vertex for mixing
    Int_t fNXEBins;//Number of bins for XE distribution (not finished)
    Int_t fNIDs;//Number of trigger IDs (such as CINT1B)
    Float_t fEffFitPt;//Pt cut used to divide efficiciency fits into low and high pt
    Int_t fNFitLowParam;//Number of fit parameters for low pt efficeincy fit
    Int_t fNFitHighParam;//Number of fit parameters for high pt efficeincy fit
    Int_t fMCHistos;//Toggle whether to make Monty Carlo histograms (if made and ran on data will simply be empty not cause a crash)

    TF1 *fFitLow;//Low pt efficiency fit function
    TF1 *fFitHigh;//High pt efficiency fit function
    Float_t *fFitLowParam; //fit parameters for low pt efficiency fit function
    Float_t *fFitHighParam;//fit parameters for high pt efficiency fit function

    Float_t *fPtTrigArray;//Array for trigger pt bins
    Float_t *fPtAssocArray;//Array for associated pt bins (2-particle correlations)
    //2 arrays to allow overlaping bins
    Float_t *fPtAssoc3Array1;//lower bin edge for 3-particle correlations associated pt bins 
    Float_t *fPtAssoc3Array2; //upper bin eged for 3-particle correlations associated pt bins
    //2 arrays to allow oeverlaping bins 
    Int_t *fCentArrayMin;//lower bin edge for centrality bins (based on SPD tracklet multiplicity)
    Int_t *fCentArrayMax;//upper bin edge for centraily bins
    Float_t *fXEArray;//XE bin array
    char *fTrigIDArray;//array of trigger ids (such as CINT1B)
    
    Float_t fVertexArray[(kNumberOfVertexBins+1)];//array for dividing z vertex into bins for event mixing
  
  //Histograms
    TH1F *fHistPt[kNumberOfCentBins][2];//Pt Distribution of tracks (no efficency correction)
    TH1F *fHistPtEff[kNumberOfCentBins][2];//Pt Distribution of tracks (efficnecy corrected if set)
    TH1F *fHistPtTrig[kNumberOfPtBins][kNumberOfCentBins][2];//Pt Distributions of tracks from triggered events
    TH1F *fHistMult[2];//Multipuliciy Distributions
    TH1F *fHistMultTrig[kNumberOfPtBins][2];//Multiplicity Distributions of triggered events

    TH2F *fHistPhi[kNumberOfCentBins][2];//Phi Distributions of tracks
    TH2F *fHistPhiTrig[kNumberOfPtBins][kNumberOfCentBins][2];//Phi Distributions fo tracks in triggered events
    TH2F *fHistDeltaPhi[kNumberOfPtBins][kNumberOfCentBins][3][2];//DeltaPhi 2-particle distributions
    TH2F *fHistDeltaPhiMix[kNumberOfPtBins][kNumberOfCentBins][3][2];//DeltaPhi 2-particle distributions from mixed events

    TH2F *fHistPhiPt[kNumberOfCentBins][2];//Pt weighed phi distributions
    TH2F *fHistPhiTrigPt[kNumberOfPtBins][kNumberOfCentBins][2];//Pt weighted phi distributions from triggered events
    TH2F *fHistDeltaPhiPt[kNumberOfPtBins][kNumberOfCentBins][3][2];//DeltaPhi 2-partricle pt weighted distributions
    TH2F *fHistDeltaPhiMixPt[kNumberOfPtBins][kNumberOfCentBins][3][2];//DeltaPhi 2-particle pt weigthed distributions from mixed events

    TH2F *fHistEta[kNumberOfCentBins][2];//Eta distributions
    TH2F *fHistEtaTrig[kNumberOfPtBins][kNumberOfCentBins][2];//Eta distributions from mixed events

    TH2F *fHistEtaPt[kNumberOfCentBins][2];//Pt-weighted Eta distributions
    TH2F *fHistEtaTrigPt[kNumberOfPtBins][kNumberOfCentBins][2];//Pt-weighted Eta diestibutions from mixed events

    TH2F *fHistDeltaEtaN[kNumberOfPtBins][kNumberOfCentBins][3][2];//Near-side 2-particle DeltaEta distributions
    TH2F *fHistDeltaEtaNMix[kNumberOfPtBins][kNumberOfCentBins][3][2];//Near-side 2-particle DeltaEta distributions from mxied events

    TH2F *fHistDeltaEtaNPt[kNumberOfPtBins][kNumberOfCentBins][3][2];//Near-side pt weighted 2-particle DeltaEta distributions
    TH2F *fHistDeltaEtaNMixPt[kNumberOfPtBins][kNumberOfCentBins][3][2];//Near-side pt-weighted 2-particle DeltaEta distributions from mixed events

    TH2F *fHistDeltaEtaA[kNumberOfPtBins][kNumberOfCentBins][3][2];//Away-side DeltaEta distributions
    TH2F *fHistDeltaEtaAMix[kNumberOfPtBins][kNumberOfCentBins][3][2];//Away-side DeltaEta distributions from mixed events

    TH2F *fHistDeltaEtaAPt[kNumberOfPtBins][kNumberOfCentBins][3][2];//Away-side pt-weighted DeltaEta distributions from mixed events
    TH2F *fHistDeltaEtaAMixPt[kNumberOfPtBins][kNumberOfCentBins][3][2];//Away-side pt-weigthed DeltaEta distributions from mixed events

    TH1F *fHistNEvents[kNumberOfCentBins][2];//Number of events count
    TH1F *fHistNTrigger[kNumberOfCentBins][2];//Number of triggers count
    TH1F *fHistNTriggerPt[kNumberOfCentBins][2];//Pt-weighted number of triggers count
    TH1F *fHistNMix[kNumberOfCentBins][2];//Number of mixed events count

    TH3F *fHistPhiEta[kNumberOfCentBins][2];//Phi-Eta distributions
    TH3F *fHistPhiEtaTrig[kNumberOfPtBins][kNumberOfCentBins][2];//Phi-Eta distributions in triggered events
    TH3F *fHistDeltaPhiEta[kNumberOfPtBins][kNumberOfCentBins][2];//DeltaPhi-DeltaEta 2-particle distributions
    TH3F *fHistDeltaPhiEtaMix[kNumberOfPtBins][kNumberOfCentBins][2];//DeltaPhi-DeltaEta 2-particle distributions from mixed events

    TH3F *fHistPhiEtaPt[kNumberOfCentBins][2];//Pt-weighted Phi-Eta distributions
    TH3F *fHistPhiEtaTrigPt[kNumberOfPtBins][kNumberOfCentBins][2];//Pt-weighted Phi-Eta distributions from triggered events
    TH3F *fHistDeltaPhiEtaPt[kNumberOfPtBins][kNumberOfCentBins][2];//Pt-weighted DeltaPhi-DeltaEta 2-particle distributions
    TH3F *fHistDeltaPhiEtaMixPt[kNumberOfPtBins][kNumberOfCentBins][2];//Pt-weigthed DeltaPhi-DeltaEta 2-particle distributions from mixed events

    TH1F *fHistXEN[kNumberOfPtBins][kNumberOfCentBins][2];//XE Near-side (not finished)
    TH1F *fHistXENMix[kNumberOfPtBins][kNumberOfCentBins][2];//XE near-side from mixed events (not finished)
    TH1F *fHistXEA[kNumberOfPtBins][kNumberOfCentBins][2];//XE away-side (not finished)
    TH1F *fHistXEAMix[kNumberOfPtBins][kNumberOfCentBins][2];//XE away-side from mxied events (not finished)

  //Three-Particle Histograms
    TH2F *fHistDeltaPhiPhi[kNumberOfPtBins][kNumberOfApt3Bins][kNumberOfCentBins][4][2];//DeltaPhi-DeltaPhi 3-particle distributions
    TH2F *fHistDeltaPhiPhiMix[kNumberOfPtBins][kNumberOfApt3Bins][kNumberOfCentBins][4][2];//DeltaPhi-DeltaPhi 3-particle distributions from mixed events
    TH2F *fHistDeltaPhiPhiSS[kNumberOfPtBins][kNumberOfApt3Bins][kNumberOfCentBins][4][2];//DeltaPhi-DeltaPhi 3-particle soft-soft background term (from mixing)
    TH2F *fHistDeltaEtaEta[kNumberOfPtBins][kNumberOfApt3Bins][kNumberOfCentBins][4][2];//DeltaEta-DeltaEta 3-particle distributions
    TH2F *fHistDeltaEtaEtaMix[kNumberOfPtBins][kNumberOfApt3Bins][kNumberOfCentBins][4][2];//DeltaEta-DeltaEta 3-particle distributions from Mixed Events
    TH2F *fHistDeltaEtaEtaSS[kNumberOfPtBins][kNumberOfApt3Bins][kNumberOfCentBins][4][2];//DeltaEta-DeltaEta 3-particle distributions soft-soft background term (from mixing)


 

  //Arrays for event mixing (only need parameters of tracks passing cuts kept)
    Float_t *fMPt[kNumberOfEventsToMix][kNumberOfCentBins][kNumberOfVertexBins][2];//Pt array (mixed events)
    Float_t *fMPhi[kNumberOfEventsToMix][kNumberOfCentBins][kNumberOfVertexBins][2]; //Phi array (mixed events)
    Int_t fMixTrack[kNumberOfEventsToMix][kNumberOfCentBins][kNumberOfVertexBins][2];//Number of tracks in the event stored (mixed events)
    Float_t *fMEta[kNumberOfEventsToMix][kNumberOfCentBins][kNumberOfVertexBins][2];//Eta Array (mixed events)
    Short_t *fMCharge[kNumberOfEventsToMix][kNumberOfCentBins][kNumberOfVertexBins][2];//Charge array (mixed events)
    Float_t *fMEff[kNumberOfEventsToMix][kNumberOfCentBins][kNumberOfVertexBins][2];//Efficiency correction array (mixed events)
    Int_t fMixPointer[kNumberOfCentBins][kNumberOfVertexBins][2];//which event should be replaced with newest mixed event (mixed events)
    Int_t fMixEnd[kNumberOfCentBins][kNumberOfVertexBins][2]; //current depth of mixed event pool (mixed events)
  //additional ones to speed up the 3-particle correlations
  Short_t *fMPtAssoc3[kNumberOfEventsToMix][kNumberOfCentBins][kNumberOfVertexBins][2][10];//which 3-particle pt bin(s) the track is in, using a bit extra memory fixing the pt bins to 10, but makes the arrays much easier (mixed events)
  Short_t *fMNPtAssoc3[kNumberOfEventsToMix][kNumberOfCentBins][kNumberOfVertexBins][2];//number of 3-particle pt bins the track is in (mixed events)
  

  //Arrays for Main Events
  Float_t *ftPhi;//Phi array (current event)
  Float_t *ftEta;//Eta array (current event)
  Float_t *ftPt;//Pt Array (current event)
  Short_t *ftCharge;//Charge Array (current event)
  Float_t *ftEff;//Efficneicy Array (current event)
  Int_t **ftPtAssoc3;//which 3-particle bins the track is in (current event)
  Int_t *ftNPtAssoc3;//number of 3-particle bins the track is in (current event)

  AliAnalysisTaskDiHadron(const AliAnalysisTaskDiHadron&);//not implimented
  AliAnalysisTaskDiHadron& operator=(const AliAnalysisTaskDiHadron&);//not implimnted

  ClassDef(AliAnalysisTaskDiHadron,1);
  
};
#endif
  
  
