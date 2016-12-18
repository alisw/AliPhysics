/// \class AliJJetJtAnalysis
/// \brief Class for jT analysis in jets
///
/// This is a longer description of the class. This longer description is
/// formatted using the Markdown syntax (see below) and can span on multiple
/// lines.
/// 
/// \author Tomas Snellman, tsnellma@cern.ch
/// \author Beomkyu Kim
/// \author Dongjo Kim
/// \date Nov 11, 2016


#ifndef ALIJJETJTANALYSIS_H
#define ALIJJETJTANALYSIS_H

#include <vector>
#include <TObjArray.h>
#include <TVector.h>
#include "AliJCard.h"
#include "AliJJet.h"
#include "AliJHistManager.h"
#include "AliJJetAnalysis.h"
#include <TRandom3.h>  //FK//

class AliJEfficiency;
class TString;

class AliJJetJtAnalysis{
  public:

    enum { kJUndefined = -1 , kJRecoTrack, kJMCParticle };

    AliJJetJtAnalysis();
    AliJJetJtAnalysis( AliJCard * card );
    AliJJetJtAnalysis(const AliJJetJtAnalysis& ap);
    AliJJetJtAnalysis& operator = (const AliJJetJtAnalysis& ap);
    ~AliJJetJtAnalysis();

    //type = charged jet or full jet, js = cone size 

    void FillHistosJets();

    void AddJets(TObjArray * jets, int type ){ 
      if( !jets ) {
        //return;
      }
      fJetListOfList.Add( (TObject*)jets ); 
      fTrackOrMCParticle.push_back(type);
      //if( !jets ) return;
      for( int i=0;i<jets->GetEntriesFast();i++ ){
        //((AliJJet*)jets->At(i))->ReSum();
      }
    } // TODO clean before event

    void AddMCJets(TObjArray * jets ){ 
      if( !jets ) {
        //return;
      }
      fMCJetListOfList.Add( (TObject*)jets ); 
      //if( !jets ) return;
      for( int i=0;i<jets->GetEntriesFast();i++ ){
        //((AliJJet*)jets->At(i))->ReSum();
      }
    } // TODO clean before event
    void SetJTracks(TClonesArray *tracks ){fTracks = tracks ;}
    void SetMCJTracks(TClonesArray *tracks ){fMCTracks = tracks ;}


    int GetNJets(){ return GetJetList()->GetEntriesFast(); }
    TObjArray* GetJetList(){ return fJetList; }
    //Double_t GetJetEtaRange(){ return fJetEtaRange; }
    //void SetJetEtaRange(double eta){ fJetEtaRange=eta; }
    void SetJetList(TObjArray* jetlist){ fJetList=jetlist; }
    void SetMCJetList(TObjArray* jetlist){ fMCJetList=jetlist; }
    void SetInputList(TObjArray * ilist){ fInputList = ilist;}
    void SetTrackOrMCParticle( UInt_t i, int v ){ fTrackOrMCParticle[i] = v; }
    int  GetTrackOrMCParticle( UInt_t i ){ return fTrackOrMCParticle.at( i ); }
    //void SetTrackJetMap(std::vector<int> * v){ fTrackJetMap=v;}
    //void SetJetPtBin( TVector * b){ fJetPtBins=b; }
    void SetCard(AliJCard * card) {fCard = card;}

    AliJJet & Jet( int i ){ return *(AliJJet*) fJetList->At(i); }
    AliJJet & MCJet( int i ){ return *(AliJJet*) fMCJetList->At(i); }

    void UserCreateOutputObjects();
    void UserExec();

    void ClearBeforeEvent();
    void CreateMCHistograms();
    void FillJtHistogram( TObjArray *Jets, int iContainer , int mc);
    void FillRandomBackground(double jetpT, double jetE, TObjArray *Jets , int iContainer, int mc);

    void FillCorrelation( TObjArray *Jets, TObjArray *mcJets, int iContainer, int iContainerParticle);
    void FillPythia(TObjArray *Jets, int iContainer); 
    int FindPythiaJet(int iContainer, int itrack);

    void FillBgJtWithSmallerR(const TClonesArray &Jets, 
        double nR, int iHist);

    void FillBgJtWithDiffAxes (int iao, int ia, int iHist);

    void WriteHistograms();

    double getDiffR(double phi1, double phi2, double eta1, double eta2);
    int GetBin(TVector *array, double val){

      int iBin=-1;
      for(int i=1; i< array->GetNoElements(); i++){
        if((*array)[i] <= val && val<(*array)[i+1]){
          iBin=i-1;
          break;
        }
      }

      return iBin;
    }
    void SetJetFinderName(vector<TString> JetFinderName){ fJetFinderName = JetFinderName; }
    void SetMCJetFinderName(vector<TString> JetFinderName){ fMCJetFinderName = JetFinderName; }

    // Need for event loop
    void SetCentralityBin( int cbin) { cBin = cbin;}
    void SetCentrality( float cent) { fcent = cent;}
    void SetZVertex( float zvtx) { zVert = zvtx;}
    void SetZVertexBin( int zbin) { zBin = zbin;}
    void SetNrandom(int nRndm) { Nrandom = nRndm;}
    void SetMoveJet(int move) { moveJet = move;}
    void SetMC(int mc) {fDoMC = mc;};
    void SetNumberOfJetFinders( int njfinder ) { nJetContainer = njfinder;}
    AliJEfficiency* GetAliJEfficiency() { return fEfficiency;}




  private:
    TObjArray * fInputList; // comment needed
    TObjArray * fJetList; ///< List of jets
    TObjArray * fMCJetList; ///< List of MC jets
    TObjArray fJetListOfList; //!<! List of jet finders
    TObjArray fMCJetListOfList; //!<! List of MC jet finders 
    vector<int>               fTrackOrMCParticle; ///< Keep track of 
    vector<TClonesArray>      fJetBgListOfList;
    TClonesArray     fpythiaJets;
    double fJetEtaCut;
    TRandom3 *frandom; // comment me

    TVector  *fJetTriggPtBorders; ///< Jet pT bin borders
    TVector  *fJetConstPtLowLimits; ///<  Comment needed
    TVector  *fJetAssocPtBorders; ///< Jet constituent pT bin borders
    TVector  *fDeltaRBorders;
    TVector *fJetLeadPtBorders; ///< Leading track pT bin borders
    TVector *fJetMultBorders; ///< Jet multiplicity bin borders
    int nJetContainer; ///< Number of jet finders

    AliJCard * fCard; // comment needed
    AliJJetAnalysis *fJJetAnalysis;
    vector<TString> fJetFinderName;
    vector<TString> fMCJetFinderName;
    vector<double>  fConeSizes;
    // Need for events
    AliJEfficiency *fEfficiency;
    int cBin;
    float fcent;
    int zBin;
    float zVert;
    TClonesArray *fTracks;
    TClonesArray *fMCTracks;
    int iConst;

    TVector *fTrackJt; ///< Store track jT values in jet cone
    TVector *fConstJt; ///< Store constituent jT values
    TVector *fTrackPt; ///< Store track pT values in jet cone
    TVector *fConstPt; ///< Store constituent jT values
    TVector *fConstLabels; 
    TVector *fJetPt; ///< Store jet pT values
    TVector *fJetPt2; ///< Comment needed
    TVector *fTrackFound; ///< Keep track of which tracks were matched with MC tracks
    TVector *fConstFound; ///< Keep track of which constituents were matched with MC tracks
    TVector *fBin2;
    TVector *fBin3;
    TVector *fpta;
    TVector *fptt;
    int Nrandom;
    int moveJet;
    int fDoMC;
    //Histograms
    AliJHistManager * fHMG;
    AliJHistManager * fHMGMC;

    AliJBin fJetFinderBin;
    AliJBin fJetTriggerBin; 
    AliJBin fTrkPtBin; 
    AliJBin fTrkLimPtBin; 
    AliJBin fJetLeadPtBin;
    AliJBin fJetMultBin;
    AliJBin fdRBin;
    AliJBin fiHist;
    AliJBin fJetFinderBinMC; 
    AliJBin fJetTriggerBinMC; 
    AliJBin fTrkPtBinMC; 
    AliJBin fTrkLimPtBinMC; 
    AliJBin fJetLeadPtBinMC;
    AliJBin fJetMultBinMC;
    AliJBin fdRBinMC;
    AliJBin fiHistMC;
    AliJBin fiHist2MC;
    AliJTH1D fhNumber;
    AliJTH1D fhKNumber;
    AliJTH1D fhJetPt; /// Jet pT distribution
    AliJTH1D fhJetPtBin; /// Jet pT distribution in jet pT bins
    AliJTH1D fhJetPtMultiplicityCutBin; /// Jet pT distribution in jet multiplicity bins
    AliJTH1D fhJetPtTrackCutBin; /// Jet pT distribution in leading track pT bins
    AliJTH1D fhJetPtWeight ; /// Jet pT distribution with \f$ \frac{1}{p_T} \f$ weight
    AliJTH1D fhJetPtWeightBin; /// Jet pT distribution with \f$ \frac{1}{p_T} \f$ weight in jet pT bins
    AliJTH1D fhJetMultiplicityBin; /// Jet multiplicity distribution from jet->GetConstituents()
    AliJTH1D fhJetMultiplicityBin2; /// Jet multiplicity based on number of tracks included in the jT distribution
    AliJTH1D fhJetConeMultiplicityBin; /// Number of tracks in constant cone around jet axis
    AliJTH1D fhZ; /// Z distribution in jets
    AliJTH1D fhZBin; /// Z distribution in jets in jet pT bins
    AliJTH1D fhRBin; /// R distribution in jets in jet pT bins \f$ R = \sqrt{\left(\Delta \phi\right)^2 + \left(\Delta \eta\right)^2}
    AliJTH1D fhJt; /// Total jT distribution in jets
    AliJTH1D fhJtBin; /// jT distribution in jets in jet pT bins
    AliJTH1D fhJtWeightBin; /// jT distribution in jets in jet pT bins weighted with \f$ \frac{1}{j_T} \f$
    AliJTH1D fhJtWeightBinTest; 
    AliJTH1D fhJtWeightBinTest2;
    AliJTH1D fhLogJtWeightBin; /// log(jT) distribution in jets in jet pT bins weighted with \f$ \frac{1}{j_T} \f$
    AliJTH1D fhLogJtWeight2Bin; /// log(jT) distribution in jets in jet pT bins weighted with \f$ \frac{1}{j_T^2} \f$
    AliJTH1D fhJtWithPtCutWeightBinBin; /// jT distribution in jets in jet pT and track pT bins weighted with \f$ \frac{1}{j_T} \f$
    AliJTH1D fhLogJtWithPtCutWeightBinBin; /// log(jT) distribution in jets in jet pT and track pT bins weighted with \f$ \frac{1}{j_T} \f$
    AliJTH1D fhLogJtWithPtCutWeight2BinBin; /// Jt distribution in jets in jet pT and track pT bins weighted with \f$ \frac{1}{j_T^2} \f$
    AliJTH1D fhJtBinLimBin; /// Comment needed
    AliJTH1D fhJtWeightBinLimBin; /// Comment needed 
    AliJTH1D fhLogJtWeightBinLimBin;
    AliJTH1D fhLogJtWeight2BinLimBin;

    //Histograms for jt in cone
    AliJTH1D fhJetConeTrkPt; /// Distribution of track pT in constant cone around jet axis
    AliJTH1D fhJetConeTrkPtBin; /// Distribution of track pT in constant cone around jet axis in jet pT bins
    AliJTH1D fhJetConeTrkPtWeightBin; /// Distribution of track pT in constant cone around jet axis in jet pT bins with \f$\frac{1}{p_T} \f$ weight 
    AliJTH1D fhJetConeZ; ///  Distribution of track Z in constant cone around jet axis \f$ z= \frac{p_{T,track}}{p_{T,jet}} \f$
    AliJTH1D fhJetConeZBin; ///  Distribution of track Z in constant cone around jet axis \f$ z= \frac{p_{T,track}}{p_{T,jet}} \f$ in jet pT bins
    AliJTH1D fhJetConeJt; /// Total jT distribution for tracks inside constant cone around jet axis
    AliJTH1D fhJetConeJtBin; /// jT distribution for tracks inside constant cone around jet axis in jet pT bins
    AliJTH1D fhJetConeJtWeightBin; /// jT distribution with \f$ \frac{1}{j_T} \f$ weight for tracks inside constant cone around jet axis in jet pT bins
    AliJTH1D fhJetConeJtWeightBinTest;
    AliJTH1D fhJetConeJtWeightBinTest2;
    AliJTH1D fhJetConeJtWeightWithTrackCutBinBin; /// jT distribution with \f$ \frac{1}{j_T} \f$ weight for tracks inside constant cone around jet axis in jet pT bins and leading track pT bins
    AliJTH1D fhJetConeJtWeightWithMultiplicityCutBinBin; /// jT distribution with \f$ \frac{1}{j_T} \f$ weight for tracks inside constant cone around jet axis in jet pT and jet multiplicity bins
    AliJTH1D fhJetConeLogJtWeightBin; /// log(jT) distribution with \f$ \frac{1}{j_T} \f$ weight for tracks inside constant cone around jet axis in jet pT bins
    AliJTH1D fhJetConeLogJtWeight2Bin; /// log(jT) distribution with \f$ \frac{1}{j_T^2} \f$ weight for tracks inside constant cone around jet axis in jet pT bins 
    AliJTH1D fhJetConeJtWithPtCutWeightBinBin; /// jT distribution with \f$ \frac{1}{j_T} \f$ weight for tracks inside constant cone around jet axis in jet pT and track pT bins
    AliJTH1D fhJetConeLogJtWithPtCutWeightBinBin; /// log(jT) distribution with \f$ \frac{1}{j_T} \f$ weight for tracks inside constant cone around jet axis in jet pT and track pT bins
    AliJTH1D fhJetConeLogJtWithPtCutWeight2BinBin; /// log(jT) distribution with \f$ \frac{1}{j_T^2} \f$ weight for tracks inside constant cone around jet axis in jet pT and track pT bins

    //Unfolding Background histograms
    AliJTH1D fhJetConeJtUnfBg; /// jT distribution of tracks inside jet cone that have no corresponding track in particle level set
    AliJTH1D fhJetConeJtBinUnfBg; /// jT distribution of tracks inside jet cone that have no corresponding track in particle level set in jet pT bins
    AliJTH1D fhJetConeJtWeightBinUnfBg; /// jT distribution with \f$\frac{1}{j_T}\f$ weight of tracks inside jet cone that have no corresponding track in particle level set in jet pT bins
    AliJTH1D fhJetConeLogJtWeightBinUnfBg; /// log(jT) distribution with \f$\frac{1}{j_T}\f$ weight of tracks inside jet cone that have no corresponding track in particle level set in jet pT bins
    AliJTH1D fhJetConeLogJtWeight2BinUnfBg; /// log(jT) distribution with \f$\frac{1}{j_T^2}\f$ weight of tracks inside jet cone that have no corresponding track in particle level set in jet pT bins

    AliJTH1D fhJetBgPt; ///
    AliJTH1D fhJetBgPtBin;
    AliJTH1D fhBgRBin; /// R distribution in background in jet pT bins \f$ R = \sqrt{\left(\Delta \phi\right)^2 + \left(\Delta \eta\right)^2}
    AliJTH1D fhBgZ; ///Background z distribution, where \f$ z = p_{track}/p_{jet} \f$. Jet refers to the jet where inclusive jT was calculated
    AliJTH1D fhBgZBin; ///Background z distribution in jet pT bins, where \f$ z = p_{track}/p_{jet} \f$. Jet refers to the jet where inclusive jT was calculated
    AliJTH1D fhBgJt; /// Background jT using perpendicular cone. Normalize with BgTrkNumberBin[iR][iJ]->GetEntries()
    AliJTH1D fhBgJtBin; ///Background jT using perpendicular cone in Jet pT bins. Normalize with BgTrkNumberBin[iR][iJ]->GetEntries()
    AliJTH1D fhBgJtWeightBin; ///Background jT weighted by \f$ \frac{1}{j_T} \f$ using perpendicular cone in Jet pT bins. Normalize with BgTrkNumberBin[iR][iJ]->GetEntries()
    AliJTH1D fhBgLogJtWeightBin; ///Background log(jT) weighted by \f$ \frac{1}{j_T} \f$ using perpendicular cone in Jet pT bins. Normalize with BgTrkNumberBin[iR][iJ]->GetEntries()
    AliJTH1D fhBgLogJtWeight2Bin; ///Background log(jT) weighted by \f$ \frac{1}{j_T^2} \f$ using perpendicular cone in Jet pT bins. Normalize with BgTrkNumberBin[iR][iJ]->GetEntries(). Difference ??? //FIXME
    AliJTH1D fhBgJtWithPtCutWeightBinBin; ///Background jT weighted by \f$ \frac{1}{j_T} \f$ using perpendicular cone in Jet pT bins and in track pT bins. Normalize with BgTrkNumberBin[iR][iJ]->GetEntries()
    AliJTH1D fhBgLogJtWithPtCutWeightBinBin; ///Background log(jT) weighted by \f$ \frac{1}{j_T} \f$ using perpendicular cone in Jet pT bins and in track pT bins. Normalize with BgTrkNumberBin[iR][iJ]->GetEntries()
    AliJTH1D fhBgLogJtWithPtCutWeight2BinBin; ///Background log(jT) weighted by \f$ \frac{1}{j_T^2} \f$ using perpendicular cone in Jet pT bins and in track pT bins. Normalize with BgTrkNumberBin[iR][iJ]->GetEntries()
    AliJTH1D fhBgJtWithPtCutWeightBinBinSmallerR;
    AliJTH1D fhBgLogJtWithPtCutWeightBinBinSmallerR;
    AliJTH1D fhBgLogJtWithPtCutWeight2BinBinSmallerR;
    AliJTH1D fhBgJtWithPtCutWeightBinBinDiffR;
    AliJTH1D fhBgLogJtWithPtCutWeightBinBinDiffR;
    AliJTH1D fhBgLogJtWithPtCutWeight2BinBinDiffR;
    AliJTH1D fhBgJtBinLimBin;
    AliJTH1D fhBgJtWeightBinLimBin;
    AliJTH1D fhBgLogJtWeightBinLimBin;
    AliJTH1D fhBgLogJtWeight2BinLimBin;
    AliJTH1D fhTrkPt;
    AliJTH1D fhTrkPtBin;
    AliJTH1D fhTrkPtWeightBin;
    AliJTH1D fhLeadingTrkPtBin;
    AliJTH1D fhBgTrkPt;
    AliJTH1D fhBgTrkPtBin;
    AliJTH1D fhBgTrkPtWeightBin;
    AliJTH1D fhBgTrkNumber;
    AliJTH1D fhBgTrkNumberBin;


    //Randomized background histograms
    AliJTH1D fhBgRndmTrkPt;
    AliJTH1D fhBgRndmZ;
    AliJTH1D fhBgRndmRBin; /// R distribution in random background in jet pT bins \f$ R = \sqrt{\left(\Delta \phi\right)^2 + \left(\Delta \eta\right)^2}
    AliJTH1D fhBgRndmJt;
    AliJTH1D fhBgRndmLogJt;
    AliJTH1D fhBgRndmJtWithPtCutWeightBin;
    AliJTH1D fhBgRndmLogJtWithPtCutWeight2Bin;
    AliJTH1D fhBgRndmJtWithPtCutWeightBinBin;
    AliJTH1D fhBgRndmLogJtWithPtCutWeight2BinBin;
    AliJTH1D fhBgRndmTrkNumber;

    AliJTH1D fhdeltaE;
    AliJTH1D fhdeltaN;
    AliJTH1D fhFullJetEChJetBin;
    AliJTH1D fhFullChdRChJetBin;
    AliJTH2D fh2DFullEvsChEdN0;
    AliJTH2D fh2DFullEvsChEdNnot0;
    AliJTH2D fhJetEtaPhi;
    AliJTH2D fhTrackEtaPhi;
    AliJTH1D fhJetdPt;

    //Jet Correlation Histograms
    AliJTH2D fhTrackJtCorrBin;	
    //AliJTH3D fhTrackJtCorr2D; //TODO
    AliJTH2D fhTrackJtCorrBinTest;	
    AliJTH2D fhTrackJtCorrBinTest2;
    AliJTH2D fhTrackPtCorr;	
    AliJTH2D fhConstPtCorr;	
    AliJTH2D fhJetPtCorr;	
    AliJTH2D fhJetPtCorr2;	
    AliJTH2D fhJetPtCorr3;	
    AliJTH2D fhJetPtCorrCoarse;	
    AliJTH1D fhJetdR;
    AliJTH1D fhTrackMatchSuccess;
    AliJTH2D fhConstJtCorrBin;	
    AliJTH1D fhConstMatchSuccess;

    //PYthia correlation histograms
    AliJTH2D fhTrackJtCorrBinPythia; /// Track jT correlation in jet pT bins between Pythia jets and measured jets
    AliJTH2D fhTrackPtCorrPythia;	/// Track pT Correlation between tracks in Pythia jets and measured jets
    AliJTH2D fhJetPtCorrPythia;	 /// Jet pT Correlation between Pythia jet candidates and measured jets
    AliJTH2D fhJetPtCorr2Pythia;	
    AliJTH2D fhJetPtCorrPythiaCoarse;	
    AliJTH1D fhJetdRPythia; /// Distance R between Pythia jet candidate and reconstructed jet \f$ R = \sqrt{\left(\Delta \phi\right)^2 + \left(\Delta \eta\right)^2} \f$
    AliJTH1D fhTrackMatchSuccessPythia;
    AliJTH1D fhJetdPtPythia; /// Difference between Pythia jet pT and reconstructed jet pT

    //Pythia jets histograms
    AliJTH1D fhJetPtPythia; /// pT distribution of Pythia jet candidates
    AliJTH1D fhJetPtBinPythia; /// pT distribution of Pythia jet candidates in jet pT bins
    AliJTH1D fhJetMultiplicityBinPythia; /// Multiplicity distribution of pythia jet fragments in jet pT bins
    AliJTH1D fhRBinPythia; /// R distribution in pythia \f$ R = \sqrt{\left(\Delta \phi\right)^2 + \left(\Delta \eta\right)^2} \f$
    AliJTH1D fhZPythia; /// z distribution in pythia \f$ z = \frac{p_{track}}{p_{grandparent}} \f$
    AliJTH1D fhZBinPythia; /// z distribution in jet pT bins in pythia \f$ z = \frac{p_{track}}{p_{grandparent}} \f$
    AliJTH1D fhZBinPythiaRCut; /// z distribution in jet pT bins in pythia \f$ z = \frac{p_{track}}{p_{grandparent}} \f$ where \f$ R < 0.4 \f$
    AliJTH1D fhJtPythia; /// jT distribution in pythia
    AliJTH1D fhJtBinPythia; /// jT distribution in pythia in jet pT bins
    AliJTH1D fhJtBinPythiaRCut; /// jT distribution in pythia in jet pT bins for tracks in R = 0.4 cone
    AliJTH1D fhJtWeightBinPythia; /// jT distribution in pythia in jet pT bins weighted by \f$ \frac{1}{j_T} \f$
    AliJTH1D fhLogJtWeightBinPythia;  /// log(jT) distribution in pythia in jet pT bins weighted by \f$ \frac{1}{j_T} \f$
    AliJTH1D fhLogJtWeight2BinPythia; /// log(jT) distribution in pythia in jet pT bins weighted by \f$ \frac{1}{j_T^2} \f$
    AliJTH1D fhJtWithPtCutWeightBinBinPythia; /// jT distribution in pythia in jet pT bins and in track pT bins weighted by \f$ \frac{1}{j_T} \f$
    AliJTH1D fhLogJtWithPtCutWeightBinBinPythia; /// log(jT) distribution in pythia in jet pT bins and in track pT bins weighted by \f$ \frac{1}{j_T} \f$
    AliJTH1D fhLogJtWithPtCutWeight2BinBinPythia; /// log(jT) distribution in pythia in jet pT bins and in track pT bins weighted by \f$ \frac{1}{j_T^2} \f$

    //Pythia background
    AliJTH1D fhBgJtPythia; /// jT background for Pythia
    AliJTH1D fhBgJtBinPythia; /// jT background for Pythia in jet pT bins
    AliJTH1D fhBgJtWeightBinPythia; /// jT background for Pythia with \f$ \frac{1}{j_T} \f$ weight in jet pT bins
    AliJTH1D fhBgLogJtWeightBinPythia; /// log(jT) background for Pythia with \f$ \frac{1}{j_T} \f$ weight in jet pT bins 
    AliJTH1D fhBgLogJtWeight2BinPythia; /// log(jT) background for Pythia with \f$ \frac{1}{j_T^2} \f$ weight in jet pT bins
    AliJTH1D fhBgJtWithPtCutWeightBinBinPythia; /// jT background for Pythia with \f$ \frac{1}{j_T} \f$ weight in jet pT bins and in track pT bins
    AliJTH1D fhBgLogJtWithPtCutWeightBinBinPythia; /// log(jT) background for Pythia with \f$ \frac{1}{j_T} weight \f$ in jet pT bins and in track pT bins
    AliJTH1D fhBgLogJtWithPtCutWeight2BinBinPythia; /// log(jT) background for Pythia with \f$ \frac{1}{j_T^2} \f$ weight in jet pT bins and in track pT bins
};

#endif

