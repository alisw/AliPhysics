#ifndef AliAnalysisTaskEmcalJetHadEPpid_h
#define AliAnalysisTaskEmcalJetHadEPpid_h

// root classes
class TClonesArray;
class TH1F;
class TH2F;
class TH3F;
class THnSparse;
class TList;
class TLorentzVector;
class TGraph;

// AliROOT classes
class AliEventPoolManager;
class AliLocalRhoParameter;
class AliEMCALTrack;
class AliMagF;
class AliESDEvent;
class AliAODEvent;
class AliEMCALGeometry;
class AliEMCALRecoUtils;
class AliESDtrack;

// container classes
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;

// includes
#include <AliAnalysisTaskEmcalJet.h>
#include <AliEmcalJet.h>
#include <AliVEvent.h>
#include <AliVTrack.h>
#include <AliVCluster.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TRandom3.h>
#include <AliLog.h>

// Local Rho includes
#include "AliAnalysisTaskLocalRho.h"
#include "AliLocalRhoParameter.h"

// PID includes
#include "AliPIDResponse.h"

#include "AliAnalysisFilter.h"

class AliAnalysisTaskEmcalJetHadEPpid : public AliAnalysisTaskEmcalJet {
 public:
  AliAnalysisTaskEmcalJetHadEPpid();
  AliAnalysisTaskEmcalJetHadEPpid(const char *name);
  //virtual ~AliAnalysisTaskEmcalJetHadEPpid() {}
  virtual ~AliAnalysisTaskEmcalJetHadEPpid();

  virtual void            UserCreateOutputObjects();
  virtual THnSparse*      NewTHnSparseF(const char* name, UInt_t entries);
  virtual void            GetDimParams(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);
  virtual THnSparse*      NewTHnSparseFPID(const char* name, UInt_t entries);
  virtual void            GetDimParamsPID(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);
  // set a bun of histogram switches up
  void                    SetPlotGlobalRho(Bool_t g)            { doPlotGlobalRho = g; } // plot global rho switch
  void                    SetVariableBinning(Bool_t v)          { doVariableBinning = v; } // do variable binning switch
  void		              SetvarbinTHnSparse(Bool_t vb)         { dovarbinTHnSparse = vb; } // variable THnSparse bin switch
  void					  SetmakeQAhistos(Bool_t QAhist)        { makeQAhistos = QAhist; } // make QA histos  
  void					  SetmakeBIAShistos(Bool_t BIAShist)    { makeBIAShistos = BIAShist; } // make bias histos
  void			          SetmakeextraCORRhistos(Bool_t Xhist)  { makeextraCORRhistos = Xhist; } // make extra correlations histos

  // set data, detectors type, and PID and PID w bias switches
  void		              SetDataType(Bool_t data)		        { useAOD = data; }    // data type switch
  void					  SetcutType(TString cut)				{ fcutType = cut; }    // EMCAL / TPC acceptance cut
  void                    SetdoPID(Bool_t p)                    { doPID = p; }   // do PID switch
  void 					  SetdoPIDtrackBIAS(Bool_t PIDbias)     { doPIDtrackBIAS = PIDbias; } // do PID track bias switch

  // give comments setter
  void					  SetdoComments(Bool_t comm)			{ doComments = comm; } // give comment switch

  // getters
  TString		  GetLocalRhoName() const		{return fLocalRhoName; }

  // set names of some objects
  virtual void            SetLocalRhoName(const char *ln)       { fLocalRhoName = ln; }
  virtual void            SetTracksName(const char *tn)         { fTracksName = tn; }
  virtual void            SetJetsName(const char *jn)           { fJetsName = jn; }

  // bias and cuts - setters
  virtual void            SetAreaCut(Double_t a)                { fAreacut    = a; }
  virtual void            SetTrkBias(Double_t b)                { fTrkBias    = b; }  //require a track with pt > b in jet
  virtual void            SetClusBias(Double_t b)               { fClusBias   = b; }  //require a cluster with pt > b in jet
  virtual void            SetTrkEta(Double_t e)                 { fTrkEta   = e; }  //eta range of the associated tracks
  virtual void            SetJetPtcut(Double_t jpt)             { fJetPtcut = jpt; } // jet pt cut
  virtual void			  SetJetRad(Double_t jrad)				{ fJetRad = jrad; } // jet radius 
  virtual void 			  SetConstituentCut(Double_t constCut)   { fConstituentCut = constCut; } // constituent Cut

  // eta and phi limits of jets - setters
  virtual void            SetJetEta(Double_t emin, Double_t emax)  { fEtamin = emin; fEtamax = emax; }
  virtual void            SetJetPhi(Double_t pmin, Double_t pmax)  { fPhimin = pmin; fPhimax = pmax; }

  // event mixing - setters
  virtual void            SetEventMixing(Int_t yesno)		   { fDoEventMixing=yesno; }
  virtual void	          SetMixingTracks(Int_t tracks)		   { fMixingTracks = tracks; }

  // jet container - setters
  void SetContainerAllJets(Int_t c)         { fContainerAllJets      = c;}
  void SetContainerPIDJets(Int_t c)         { fContainerPIDJets      = c;}

protected:
  // functions 
  void					 ExecOnce();
  Bool_t		         Run();
  virtual void           Terminate(Option_t *); 
  virtual Int_t          AcceptMyJet(AliEmcalJet *jet);   // applies basic jet tests/cuts before accepting
  virtual Int_t          GetCentBin(Double_t cent) const; // centrality bin of event
  Float_t                RelativePhi(Double_t mphi,Double_t vphi) const; // relative jet track angle
  Float_t                RelativeEPJET(Double_t jetAng, Double_t EPAng) const;  // relative jet event plane angle
  virtual Int_t          GetEtaBin(Double_t eta) const;      // eta bins
  virtual Int_t          GetpTjetBin(Double_t pt) const;     // jet pt bins
  virtual Int_t          GetpTtrackBin(Double_t pt) const;   // track pt bins
  virtual Int_t          GetzVertexBin(Double_t zVtx) const; // zVertex bin
  void                   SetfHistPIDcounterLabels(TH1* fHistPID) const;  // PID counter

  // parameters of detector to cut on for event
  Double_t               fPhimin;                  // phi min
  Double_t               fPhimax;                  // phi max
  Double_t               fEtamin;                  // eta min
  Double_t               fEtamax;                  // eta max
  Double_t               fAreacut;                 // area cut
  Double_t               fTrkBias;                 // track bias
  Double_t               fClusBias;                // cluster bias
  Double_t               fTrkEta;                  // eta min/max of tracks
  Double_t	             fJetPtcut;		           // jet pt to cut on for correlations
  Double_t				 fJetRad;				   // jet radius
  Double_t				 fConstituentCut;          // jet constituent cut

  // event mixing
  Int_t			 fDoEventMixing;
  Int_t			 fMixingTracks;

  // switches for plots
  Bool_t		 doPlotGlobalRho;
  Bool_t		 doVariableBinning;
  Bool_t         dovarbinTHnSparse;
  Bool_t		 makeQAhistos;
  Bool_t		 makeBIAShistos;
  Bool_t		 makeextraCORRhistos; 

  // data type switch
  Bool_t	     useAOD;

  // Cut type (EMCAL/TPC acceptance)
  TString        fcutType;

  // switches for PID
  Bool_t		 doPID;
  Bool_t		 doPIDtrackBIAS;

  // do comment switch
  Bool_t		 doComments;

  // local rho value
  Double_t		 fLocalRhoVal;

  // object names
  TString		 fTracksName;
  TString		 fJetsName;

  // event counter
  Int_t			 event;

  // boolean functions for PID
  Bool_t		 isPItpc, isKtpc, isPtpc;
  Double_t			 nPIDtpc;

  Bool_t		 isPIits, isKits, isPits;
  Double_t			 nPIDits;  

  Bool_t		 isPItof, isKtof, isPtof;
  Double_t			 nPIDtof;

  // event pool
  TObjArray*		    CloneAndReduceTrackList(TObjArray* tracks);
  AliEventPoolManager   *fPoolMgr;  // event pool Manager object

  // PID
  AliPIDResponse	*fPIDResponse;   // PID response object
  AliTPCPIDResponse	*fTPCResponse;   // TPC pid response object

 private:
  // needed for PID, track objects
  AliESDEvent       *fESD;          // ESD object
  AliAODEvent	    *fAOD;		  // AOD object

  TH2F                  *fHistTPCdEdX;
  TH2F	                *fHistITSsignal;
//  TH2F		    *fHistTOFsignal;

  TH2F                  *fHistRhovsCent; //!
  TH2F                  *fHistNjetvsCent;//! number of jets versus Centrality
  TH2F                  *fHistJetPtvsTrackPt[6];//!
  TH2F                  *fHistRawJetPtvsTrackPt[6];//!
  TH1F                  *fHistTrackPt[6];//!
  TH1F                  *fHistEP0[6];//!
  TH1F                  *fHistEP0A[6];//!
  TH1F                  *fHistEP0C[6];//!
  TH2F                  *fHistEPAvsC[6];//!
  TH1F					*fHistJetPtcorrGlRho[6];//!
  TH2F                  *fHistJetPtvsdEP[6];//!
  TH2F                  *fHistJetPtvsdEPBias[6];//!
  TH2F                  *fHistRhovsdEP[6]; //!
  TH3F                  *fHistJetEtaPhiPt[6];//!
  TH3F                  *fHistJetEtaPhiPtBias[6];//!
  TH2F                  *fHistJetPtArea[6];//!
  TH2F                  *fHistJetPtAreaBias[6];//!
  TH2F                  *fHistJetPtNcon[6]; //!
  TH2F                  *fHistJetPtNconBias[6]; //!
  TH2F                  *fHistJetPtNconCh[6]; //!
  TH2F                  *fHistJetPtNconBiasCh[6]; //!
  TH2F                  *fHistJetPtNconEm[6]; //!
  TH2F                  *fHistJetPtNconBiasEm[6]; //!
  TH1F		            *fHistJetHaddPhiINcent[6];
  TH1F			        *fHistJetHaddPhiOUTcent[6];
  TH1F			        *fHistJetHaddPhiMIDcent[6];

  TH1                   *fHistCentrality;
  TH1                   *fHistZvtx;
  TH1                   *fHistMult;
  TH1		        	*fHistJetPhi;
  TH1		            *fHistTrackPhi;
  TH1		            *fHistJetHaddPhiIN;
  TH1			        *fHistJetHaddPhiOUT;
  TH1			        *fHistJetHaddPhiMID;
  TH1					*fHistJetHaddPhiBias;
  TH1					*fHistJetHaddPhiINBias;
  TH1					*fHistJetHaddPhiOUTBias;
  TH1					*fHistJetHaddPhiMIDBias;

  TH1                   *fHistMEdPHI; // phi distrubtion of mixed events
  TH1					*fHistTrackPtallcent;

  TH2                   *fHistJetEtaPhi;  
  TH2                   *fHistTrackEtaPhi[4][7];
  TH1		        	*fHistJetHadbindPhi[9]; 
  TH1					*fHistJetHadbindPhiIN[9]; 
  TH1					*fHistJetHadbindPhiMID[9]; 
  TH1		       		*fHistJetHadbindPhiOUT[9]; 
  TH2                   *fHistJetHEtaPhi;

  TH1                   *fHistJetPt[2];
  TH1                   *fHistJetPtBias[2];
  TH1                   *fHistJetPtTT[2];
  TH2                   *fHistAreavsRawPt[2];
  TH2                   *fHistJetH[2][5][3];
  TH2                   *fHistJetHBias[2][5][3];
  TH2                   *fHistJetHTT[2][5][3];
  TH1F					*fHistJetHdPHI[11];
  TH2F					*fHistJetHdETAdPHI[11];
  TH2F                  *fHistSEphieta; // single events phi-eta distributions
  TH2F                  *fHistMEphieta; // mixed events phi-eta distributions
  TH1F					*fHistJetHaddPHI;

  // PID status histo's
  TH1					*fHistPID;

  // THn Sparse's
  THnSparse             *fhnPID;          // PID sparse
  THnSparse             *fhnMixedEvents;  // mixed events matrix
  THnSparse             *fhnJH;           // jet hadron events matrix

  // container objects
  AliJetContainer            *fJetsCont;                   //!Jets
  AliParticleContainer       *fTracksCont;                 //!Tracks
  AliClusterContainer        *fCaloClustersCont;           //!Clusters

  // container specifier
  Int_t					fContainerAllJets;  // number of container with all full jets
  Int_t					fContainerPIDJets;  // number of container with full jets meeting Pt cut (for PID)

// ***********************************************************
   
  //Declare it private to avoid compilation warning
  AliAnalysisTaskEmcalJetHadEPpid(const AliAnalysisTaskEmcalJetHadEPpid & g) ; // cpy ctor

  AliAnalysisTaskEmcalJetHadEPpid& operator=(const AliAnalysisTaskEmcalJetHadEPpid&); // not implemented
  ClassDef(AliAnalysisTaskEmcalJetHadEPpid, 4); // Emcal jet hadron PID - Event plane dependence
};
#endif
