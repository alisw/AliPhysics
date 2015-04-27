#ifndef ALIANALYSISTASKJETJTJT_H
#define ALIANALYSISTASKJETJTJT_H

class TH1;
class TH2;
class TH3;
class TGraphErrors;
class TProfile;
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;
class AliAnalysisUtils;
class JTJTEfficiency;


#include "AliAnalysisTaskEmcalJet.h"

using std::cout;
using std::endl;

class AliAnalysisTaskJetJTJT : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskJetJTJT();
  AliAnalysisTaskJetJTJT(const char *name);
  virtual ~AliAnalysisTaskJetJTJT();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

  void SetTrackArrayName( char *c ) { fTrackArrayName = c; }

  void setCentBinBorders( int n, Double_t *c);
  void setTriggPtBorders( int n, Double_t *c);
  void setAssocPtBorders( int n, Double_t *c);
  void setEffMode(int m) {fEffMode = m; };
  void setDebug(int n) {debug = n; }
  void setRunPeriod(const char *period) {runPeriod = period; cout << "RunPeriod set to " << runPeriod << endl;}
  //AliJetContainer            *AddJetContainer(const char *n, TString defaultCutType = "", Float_t jetRadius = 0.4);


 protected:
  void                        ExecOnce();
  Bool_t                      FillHistograms()   ;
  Bool_t                      Run()              ;
  void                        CheckClusTrackMatching();
  Double_t	     	      getJt(AliVTrack *track, AliEmcalJet *jet, int reverse);
  Double_t	     	      getJt(AliVParticle *track, AliEmcalJet *jet, int reverse);
  Double_t                    getDiffR(double phi1, double phi2, double eta1, double eta2);


  // General histograms
  TH1                       **fHistTracksPt;            //!Track pt spectrum
  TH1                       **fHistTracksJt;            //!Track jt spectrum
  TH1			    **fHistTracksEta;           //!Track eta spectrum
  TH1                       **fHistClustersPt;          //!Cluster pt spectrum
  TH1                       **fHistLeadingJetPt;        //!Leading jet pt spectrum
  TH1                       ***fHistJetsPt;          	//!Jet pt spectrum
  TH1                       ***fHistJetsCorrPt;        	//!Rho corrected Jet pt spectrum
  TProfile                  **fHistJetsCorrPtVsNonCorr; //!Corrected versus raw jet pt
  TH1			    ***fHistBackgroundDone;	//!Background test

  //Jt histograms
  TH1                       ****fHistJTPta;			//!Jet Jt spectrum
  TH1                       ****fHistLogJTPta;			//!Logarithmic Jet Jt spectrum
  TH1                       ****fHistJTPta_all;			//!All particles Jt spectrum
  TH1                       ****fHistJTBg;			//!Jt background
  TH1                       ****fHistLogJTBg;			//!Logarithmic Jt background

  TProfile		    ***fHistPtaVsJt;             //!Associated pT vs. Jt in jet
  TProfile		    ***fHistBgPtaVsJt;           //!Associated pT vs. Jt in background

  //non-invariant
  TH1                       ****fHistJTPtaNonInv;			//!Jet Jt spectrum
  TH1                       ****fHistLogJTPtaNonInv;			//!Logarithmic Jet Jt spectrum
  TH1                       ****fHistJTPta_allNonInv;			//!All particles Jt spectrum
  TH1                       ****fHistJTBgNonInv;			//!Jt background
  TH1                       ****fHistLogJTBgNonInv;			//!Logarithmic Jt background


  //Background statistics

  TH1			    ***fHistBgMulti;		//!Multiplicity in background cone
  TH1			    ***fHistBgPt;		//!Background pt distribution

  //Jet statistics
  TH1			    ***fHistJetEta;		//!Jet eta distribution
  TH1			    ***fHistJetMulti;		//!Multiplicity in jet				
  TH1			    ***fHistJetTracksPt;		//!Track pT in jet				
  TProfile		    **fhTrackingEfficiency; 	//!Tracking efficiency

  Int_t fNpttBins;
  Int_t fNptaBins;
  Int_t fEffMode;

  AliJetContainer            *fJetsCont;                   //!Jets
  //AliJetContainer            **fJetsConts;              //!Jets
  //Int_t 		     nJetsConts;
  AliParticleContainer       *fTracksCont;                 //!Tracks
  AliClusterContainer        *fCaloClustersCont;           //!Clusters  

  TH1F     *fhVertexZ;  //! vertexZ inclusive
  TH1I	   *fHistEvtSelection; //! Event selection statistics	

 private:
  //TVector *EtaGapThresholds;
  //TVector *RGapThresholds;
  //TVector *KlongBorders;
  //TVector *XeBorders;


  const AliVVertex*   fPrimaryVertex;         //! Vertex found per event


  TClonesArray *fTracks;  //! tracks array
  TString    fTrackArrayName; // track constituents array name
  TString    runPeriod; // run period name
  JTJTEfficiency *fEfficiency; //! AliJ Efficiency
  AliAnalysisUtils* fVertexHelper;	//! Vertex selection helper


  //TVector *CentBinBorders;
  //TVector *TriggPtBorders;
  //TVector *AssocPtBorders;
  Double_t CentBinBorders[10];
  Double_t TriggPtBorders[10];
  Double_t AssocPtBorders[10];
  AliAnalysisTaskJetJTJT(const AliAnalysisTaskJetJTJT&);            // not implemented
  AliAnalysisTaskJetJTJT &operator=(const AliAnalysisTaskJetJTJT&); // not implemented

  Int_t debug;

  ClassDef(AliAnalysisTaskJetJTJT, 3) // jet sample analysis task



};

class JTJTEfficiency {  // this part can occurr anywhere inside the outer accolades
#define JUNUSED(expr) do { (void)(expr); } while (0)
	public:
		enum Mode { kNotUse, kPeriod, kRunNumber, kAuto };
		enum { kJTPCOnly, kJRaa, kJGlobalTightDCA, kJGlobalDCA, kJGlobalSDD , kJHybrid, kJNTrackCuts };



		JTJTEfficiency();
		JTJTEfficiency(const JTJTEfficiency& obj);
		JTJTEfficiency& operator=(const JTJTEfficiency& obj);
		void SetMode( int i ){ fMode = i; }
		void SetDataPath(TString s ){ fDataPath=s; }
		void SetPeriodName(TString s ){ fPeriodStr=s; cout << "Eff: Run Period is set to " << fPeriodStr << endl; }
		TString GetName() const { return fName; }
		double GetCorrection( double pt, int icut, double cent ) const ;
		void SetRunNumber( Long64_t runnum ){ fRunNumber=runnum; }

		TString GetEffName() ;
		TString GetEffFullName() ;
		bool   Load();
		void   PrintOut() const {
			cout<<fInputRootName<<endl;
		}
		//void Write();

	private:
		int      fMode;             // Mode. see enum Mode
		int      fPeriod;           // Data Period index
		//AliJTrackCut fTrackCut;     // Track Cut Object. TODO:why not pointer?
		//AliJRunTable fRunTable;     // run Table. TODO:why not pointer?

		TString fDataPath;          // locaction of eff files
		TString fName;              // name of efficiency. usually empty
		TString fPeriodStr;         // DATA period
		TString fMCPeriodStr;       // MC period
		Long64_t fRunNumber;        // Runnumber
		TString fTag;               // Tags to distinguish special eff file
		TString fInputRootName;     // name of input

		TFile * fInputRoot;         // input file  
		TDirectory * fEffDir[3];    // root directory of efficiency. only second item of fEffDir with "Efficiency" is being used.
		TGraphErrors * fCorrection[20][20][20]; // Storage of Correction factor 
		TAxis * fCentBinAxis;     // Bin of Centrality. replace with AliJBin?



};

#endif
