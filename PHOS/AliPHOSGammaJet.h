#ifndef ALIPHOSGammaJet_H
#define ALIPHOSGammaJet_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.9  2005/12/20 07:31:51  schutz
 * added data members
 *
 * Revision 1.7  2005/05/28 14:19:04  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
//  Class for the analysis of gamma-jet correlations.     
//  Basically it seaches for a prompt photon in the PHOS acceptance, 
//  if so we construct a jet around the highest pt particle in the opposite 
//  side in azimuth. This jet has to fullfill several conditions to be 
//  accepted. Then the fragmentation function of this jet is constructed 

//*-- Author: Gustavo Conesa & Yves Schutz (IFIC, CERN)

// --- ROOT system ---
#include "TTask.h"
#include "TRandom.h"
#include "TArrayD.h"
#include "AliESD.h"
#include "TChain.h"

class AliPHOSFastGlobalReconstruction ;

//#include "../PYTHIA6/AliGenPythia.h"

// --- AliRoot header files ---

class AliPHOSGammaJet : public TTask {

public: 

  AliPHOSGammaJet() ; // default ctor
  AliPHOSGammaJet(const TString inputfilename) ; //ctor 
  AliPHOSGammaJet(const AliPHOSGammaJet & gj) ; // cpy ctor
  AliPHOSGammaJet & operator = (const AliPHOSGammaJet & /*rvalue*/) 
    { return *this ;} //assignement operator requested by coding convention but not needed
  virtual ~AliPHOSGammaJet() ; //virtual dtor
  virtual void   Exec(Option_t *option); 
  void List() const; 
   Double_t GetAngleMaxParam(Int_t i) const {return fAngleMaxParam.At(i) ; }
   Double_t GetEtaCut() const {return fEtaCut;}
   Double_t GetPhiEMCALCut(Int_t i) const {return fPhiEMCALCut[i];}
   TString  GetHIJINGFileName() const {return fHIJINGFileName ; }
   TString  GetHistosFileName() const {return fOutputFileName ; }
   Double_t GetInvMassMaxCut() const {return fInvMassMaxCut ; }
   Double_t GetInvMassMinCut() const {return fInvMassMinCut ; }
   Double_t GetPhiMaxCut() const {return fPhiMaxCut ; }
   Double_t GetPhiMinCut() const {return fPhiMinCut ; }
   Double_t GetPtCut() const {return fPtCut ; }
   Double_t GetNeutralPtCut() const {return fNeutralPtCut ; }
   Double_t GetChargedPtCut() const {return fChargedPtCut ; }
   Double_t GetPtJetSelectionCut() const {return fPtJetSelectionCut ; }
   Double_t GetMinDistance() const {return fMinDistance ; }
   Double_t GetJetRatioMaxCut() const {return fJetRatioMaxCut ; }
   Double_t GetJetRatioMinCut() const {return fJetRatioMinCut ; }
   Double_t GetRatioMaxCut() const {return fRatioMaxCut ; }
   Double_t GetRatioMinCut() const {return fRatioMinCut ; }
   Int_t    GetNEvent() const {return fNEvent ; }
   Int_t    GetNCones() const {return fNCone ; }
   Int_t    GetNPtThres() const {return fNPt ; }
   Float_t  GetCone() const {return fCone ; }
   Float_t  GetPtThreshold() const {return fPtThreshold ; }
   Float_t  GetCones(Int_t i) const {return fCones[i] ; }
   Float_t  GetPtThreshold(Int_t i) const {return fPtThres[i] ; }
   TString  GetConeName(Int_t i) const {return fNameCones[i] ; }
   TString  GetPtThresName(Int_t i) const {return fNamePtThres[i] ; }
   Bool_t   GetTPCCutsLikeEMCAL() const {return fTPCCutsLikeEMCAL ; }
   Bool_t   IsESDdata() const {return fESDdata ; }

   TString  GetDirName() const {return fDirName ; }
   TString  GetESDTreeName() const {return fESDTree ; }
   char *   GetDirPattern()  const {return fPattern ; }

   Bool_t   IsAnyConeOrPt() const {return fAnyConeOrPt ; }
   Bool_t   IsFastReconstruction() const {return fOptFast ; }
   Bool_t   IsHIJING() const {return fHIJING ; }
   Bool_t   IsOnlyCharged() const {return fOnlyCharged ; }

  void Plot(TString what="all", Option_t *option="") const;
  void Print(const Option_t * opt)const;

  void SetAngleMaxParam(Int_t i, Double_t par)
  {fAngleMaxParam.AddAt(par,i) ; }
  void SetAnyConeOrPt(Bool_t any){fAnyConeOrPt = any ;}
  void SetEtaCut(Double_t etacut) {fEtaCut = etacut ; }
  void SetPhiEMCALCut(Double_t phi, Int_t i){fPhiEMCALCut[i] = phi; }
  void SetFastReconstruction(Bool_t fast){fOptFast = fast ; }
  void SetHIJING(Bool_t opt){fHIJING = opt; }
  void SetHIJINGFileName(TString file){fHIJINGFileName = file ; }
  void SetNCones(Int_t n){fNCone = n ; }
  void SetNPtThresholds(Int_t n){fNPt = n ; }
  void SetCones(Int_t i, Float_t cone, TString sc)
    {fCones[i] = cone ; fNameCones[i] = sc; };
  void SetCone(Float_t cone)
    {fCone = cone; }
  void SetPtThreshold(Float_t pt){fPtThreshold = pt; };
  void SetPtThresholds(Int_t i,Float_t pt, TString spt){fPtThres[i] = pt ; 
  fNamePtThres[i] = spt; };
  void SetOnlyCharged(Bool_t opt){fOnlyCharged = opt; }
  void SetPythiaFileName(TString file){fInputFileName = file ; }
  void SetHistosFileName(TString file){fOutputFileName = file ; }
  void SetInvMassCutRange(Double_t invmassmin, Double_t invmassmax)
    {fInvMassMaxCut =invmassmax;  fInvMassMinCut =invmassmin;}	
  void SetJetRatioCutRange(Double_t ratiomin, Double_t ratiomax)
    {fJetRatioMaxCut =ratiomax;  fJetRatioMinCut = ratiomin ; }
 void SetJetTPCRatioCutRange(Double_t ratiomin, Double_t ratiomax)
    {fJetTPCRatioMaxCut =ratiomax;  fJetTPCRatioMinCut = ratiomin ; }
  void SetNEvent(Int_t n){fNEvent  = n ; }
  void SetMinDistance(Double_t min){fMinDistance  = min ; }
  void SetPhiCutRange(Double_t phimin, Double_t phimax)
  {fPhiMaxCut =phimax;  fPhiMinCut =phimin;}
  void SetPtCut(Double_t ptcut)
  {fPtCut =ptcut;}
  void SetNeutralPtCut(Double_t ptcut)
  {fNeutralPtCut =ptcut;}
  void SetChargedPtCut(Double_t ptcut)
  {fChargedPtCut =ptcut;}
  void SetPtJetSelectionCut(Double_t cut){fPtJetSelectionCut = cut; }
  void SetJetSelection(Bool_t select){ fSelect= select ; }
  void SetRatioCutRange(Double_t ratiomin, Double_t ratiomax)
  {fRatioMaxCut = ratiomax;  fRatioMinCut = ratiomin;}
  void SetTPCCutsLikeEMCAL(Bool_t b){ fTPCCutsLikeEMCAL= b ; }
  void SetESDdata(Bool_t b){ fESDdata= b ; }
 
  void  SetDirName(TString dn) { fDirName  = dn ; }
  void  SetESDTreeName(TString en) { fESDTree = en ; }
  void  SetDirPattern(char * dp)  { fPattern = dp ; }

 private:
//   void AddHIJINGToList(TList & particleList, TList & particleListCh, 
// 		       TList & particleListNe, const Int_t iEvent, 
// 		       const TLorentzVector gamma, Double_t & rot ); 
 
  void AddHIJINGToList(Int_t iEvent, TClonesArray * particleList, 
		       TClonesArray * plCh, TClonesArray * plNe, 
		       TClonesArray * plNePHOS); 


  Double_t CalculateJetRatioLimit(const Double_t ptg, const Double_t *param, 
				  const Double_t *x);

  void CreateParticleList(Int_t iEvent, TClonesArray * particleList, 
			  TClonesArray * plCh, TClonesArray * plNe, 
			  TClonesArray * plNePHOS); 
  void CreateParticleListFromESD(TClonesArray * particleList, 
				 TClonesArray * plCh, TClonesArray * plNe, 
				 TClonesArray * plNePHOS,  
				 const AliESD * esd );//, Int_t iEvent); 
  
  void FillJetHistos(TClonesArray * pl, Double_t ptg, TString conf, TString type);

  void FillJetHistosAnyConeOrPt( TClonesArray * pl, Double_t ptg, TString conf, 
				 TString type, TString cone, TString ptcut);
  Bool_t IsAngleInWindow(const Float_t angle, const Float_t e);
  Bool_t IsJetSelected(const Double_t ptg, const Double_t ptjet, 
		       const TString type);

  void MakeJet(TClonesArray * particleList, 
	       Double_t ptg, Double_t phig,
	       Double_t ptl, Double_t phil, Double_t etal, 
	       TString  type, TLorentzVector & jet); 
  void MakeJetAnyConeOrPt(TClonesArray * particleList, Double_t ptg, 
			  Double_t phig, Double_t ptl, Double_t phil, 
			  Double_t etal, TString  type); 
  void GetGammaJet(TClonesArray * pl,  Double_t &pt, 
		   Double_t &phi, Double_t &eta, Bool_t &Is)  const;

  void GetLeadingCharge(TClonesArray * pl, 
			Double_t ptg,  Double_t phig, 
			Double_t &pt, Double_t &eta, Double_t &phi) const ;
  void GetLeadingPi0   (TClonesArray * pl, 
			Double_t ptg, Double_t phig, 
			Double_t &pt, Double_t &eta, Double_t &phi)  ;

  void InitParameters();
  Double_t MakeEnergy(const Double_t energy) ;
  void MakeHistos() ;
  void MakePhoton(TLorentzVector & particle) ; 
  TVector3 MakePosition(const Double_t energy, const TVector3 pos) ;
 
  void Pi0Decay(Double_t mPi0, TLorentzVector &p0, 
		TLorentzVector &p1, TLorentzVector &p2, Double_t &angle) ;

  TChain * ReadESDfromdisk(const UInt_t eventsToRead,
			      const TString dirName, 
			      const TString esdTreeName = "esdTree",
			      const char *  pattern     = "Evt") ;
  TChain * ReadESD(const UInt_t eventsToRead = 1,
		      TString fileName  = "",
		      const TString esdTreeName = "esdTree",
		      const char *  pattern     = "Evt" ) ;
  
  Double_t SigmaE(Double_t energy) ;
  Double_t SigmaP(Double_t energy) ;
  
  void SetJet(TParticle * part, Bool_t & b, Float_t cone, Double_t eta, 
	      Double_t phi);

 private: 
  Bool_t     fAnyConeOrPt;     //  To play with the jet cone size and pt th.
  Option_t * fOptionGJ ;       //! Fill most interesting histograms 
		               // and give interesting information
  TFile *    fOutputFile ;     //! Output file
  TString    fOutputFileName;  //! Output file Name
  TString    fInputFileName;   //!
  TString    fHIJINGFileName;  //!
  Bool_t     fHIJING;          // Add HIJING event to PYTHIA event?
  Bool_t     fESDdata ;        // Read ESD?      
  Double_t   fEtaCut ;         // Eta cut
  Bool_t     fOnlyCharged ;    // Only jets of charged particles
  Double_t   fPhiEMCALCut[2] ; // Phi cut maximum
  Double_t   fPhiMaxCut ;      // Phi cut maximum
  Double_t   fPhiMinCut ;      // Phi cut minimun
  Double_t   fPtCut ;          // Min pt in PHOS
  Double_t   fNeutralPtCut ;   // Min pt detected in PHOS
  Double_t   fChargedPtCut ;   // Min pt detected in TPC
  Double_t   fInvMassMaxCut ;  // Invariant Mass cut maximum
  Double_t   fInvMassMinCut ;  // Invariant Masscut minimun
  Double_t   fMinDistance ;    // Minimal distance to resolve gamma decay.
  Double_t   fRatioMaxCut ;    // Leading particle/gamma Ratio cut maximum
  Double_t   fRatioMinCut ;    // Leading particle/gamma Ratio cut minimum
  Bool_t     fTPCCutsLikeEMCAL ; //Same jet energy ratio limits for both conf.

  //Read ESD Paramenters
  TString fDirName ; //Name of the directory where data is
  TString fESDTree ; //Name of the ESD Tree
  char *  fPattern ; //Pattern followed by  directory data

  //Jet selection parameters
  //Fixed cuts (old)
  Double_t   fJetTPCRatioMaxCut ; // Leading particle/gamma Ratio cut maximum
  Double_t   fJetTPCRatioMinCut ; // Leading particle/gamma Ratio cut minimum
  Double_t   fJetRatioMaxCut ; // Jet/gamma Ratio cut maximum
  Double_t   fJetRatioMinCut ; // Jet/gamma Ratio cut minimum

  //Cuts depending on jet pt
  Double_t fJetE1[2];    //Rec. jet energy parameters
  Double_t fJetE2[2];    //Rec. jet energy parameters
  Double_t fJetSigma1[2];//Rec. sigma of jet energy  parameters
  Double_t fJetSigma2[2];//Rec. sigma of jet energy  parameters
  Double_t fBkgMean[6];  //Background mean energy 
  Double_t fBkgRMS[6];   //Background RMS
  Double_t fJetXMin1[6]; //X Factor to set jet min limit for pp
  Double_t fJetXMin2[6]; //X Factor to set jet min limit for PbPb
  Double_t fJetXMax1[6]; //X Factor to set jet max limit for pp
  Double_t fJetXMax2[6]; //X Factor to set jet max limit for PbPb

  Int_t      fNEvent ;           // Number of events to analyze
  Int_t      fNCone ;            // Number of jet cones sizes
  Int_t      fNPt   ;            // Number of jet particle pT threshold
  Double_t   fCone  ;            // Jet cone sizes under study (!fAnyConeSize)
  Double_t   fCones[10];         // Jet cone sizes under study (fAnyConeSize)
  TString    fNameCones[10];     // String name of cone to append to histos
  Double_t   fPtThreshold;       // Jet pT threshold under study(!fAnyConeSize)
  Double_t   fPtThres[10];       // Jet pT threshold under study(fAnyConeSize)
  Double_t   fPtJetSelectionCut; // Jet pt to change to low pt jets analysis
  TString    fNamePtThres[10];   // String name of pt th to append to histos
  TObjArray * fListHistos ;      //! list of Histograms
  AliPHOSFastGlobalReconstruction * fFastRec; //Pointer to fast recons.
  Bool_t     fOptFast;    // Do we want fast Rec?
  TRandom    fRan ;       //! random number generator
  //Energy and position parameters
  Double_t   fResPara1 ;  // parameter for the energy resolution dependence  
  Double_t   fResPara2 ;  // parameter for the energy resolution dependence  
  Double_t   fResPara3 ;  // parameter for the energy resolution dependence 
  Double_t   fPosParaA ;  // parameter for the position resolution
  Double_t   fPosParaB ;  // parameter for the position resolution 
  TArrayD    fAngleMaxParam ; //Max opening angle selection parameters
  Bool_t fSelect  ;  //Select jet within limits

  ClassDef(AliPHOSGammaJet,3)
} ;
 

#endif //ALIPHOSGammaJet_H



