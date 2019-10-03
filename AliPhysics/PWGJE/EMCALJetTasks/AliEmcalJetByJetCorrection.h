#ifndef ALIEMCALJETBYJETCORRECTION_H
#define ALIEMCALJETBYJETCORRECTION_H

#include <TNamed.h>
#include <TObjArray.h>
#include <TArrayF.h>
#include <TClonesArray.h>

class TH3;
class TProfile;
class AliEmcalJet;
class TList;
class THnSparse;

class AliEmcalJetByJetCorrection : public TNamed
{
 public:

  AliEmcalJetByJetCorrection();
  AliEmcalJetByJetCorrection(const char* name);
  AliEmcalJetByJetCorrection(const AliEmcalJetByJetCorrection &jet);
  AliEmcalJetByJetCorrection& operator=(const AliEmcalJetByJetCorrection &jet);
  ~AliEmcalJetByJetCorrection() {;}

  void SetTemplate(TH3 *h)          { fh3JetPtDRTrackPt = h; }
  void SetJetPtBinWidth(Double_t w) { fBinWidthJetPt    = w; }
  void SetJetPtRange(Double_t min, Double_t max) {fJetPtMin = min; fJetPtMax = max;}
  void SetFixedTrackEfficiency(Double_t eff) { fEfficiencyFixed = eff; }
  void SetEfficiencyHist(TH1 *h)             { fhEfficiency     = h  ; }
  void SetCorrectTemplateTrackpT(Bool_t correct=kTRUE) {fCorrectpTtrack = correct;}
  void SetPoissonianNmissing(Bool_t set=kTRUE) {fNpPoisson=set;}
  void SetNMissedTracks(Double_t number) {fNMissedTracks = number;}
  void SetExternalDefinitionOfNmissed(Bool_t set=kTRUE) {fExternalNmissed=set;}
  
  Int_t              GetJetPtBin(const Double_t jetpt) const;
  Double_t           GetEfficiency(const Double_t pt) const;
  Double_t           GetMeanPtConstituents(const AliEmcalJet *jet, TClonesArray *fTracks) const;
                    
  TProfile          *GetAppliedEfficiency() const      {return fpAppliedEfficiency;}
  Bool_t             GetCorrectTemplateTrackpT() const {return fCorrectpTtrack;}
  Bool_t             GetPoissonianNmissing() const     {return fNpPoisson;}
  Bool_t             GetExternalDefinitionOfNmissed() const {return fExternalNmissed;}
  TList             *GetListOfOutput()       const     {return fListOfOutput;}
  void               Init();
  AliEmcalJet       *Eval(const AliEmcalJet *jet, TClonesArray *fTracks);
  TClonesArray      *GetAddedTrackArray()   const      {return fArrayTrackCorr;}
  
 protected:
  TH3      *fh3JetPtDRTrackPt;                 ///< 3D template
  Double_t  fBinWidthJetPt;                    ///< jet pt bin width in which to do correction
  Double_t  fJetPtMin;                         ///< minimum jet pt
  Double_t  fJetPtMax;                         ///< maximum jet pt
  TObjArray fCollTemplates;                    ///< templates (2D histos with track pT vs r)
  Bool_t    fInitialized;                      ///< status of initialization
  Double_t  fEfficiencyFixed;                  ///< fixed efficiency for all pT and all types of tracks
  TH1      *fhEfficiency;                      ///< single particle efficiency
  TH1D     *fhSmoothEfficiency;                ///< single particle efficiency smooth (see Init())
  Bool_t   fCorrectpTtrack;                    ///< if true the templates are corrected by track efficiency
  Bool_t   fNpPoisson;                         ///< draw Nmissing particle from a Poissonian with mean Nconst(1/eff-1)
  Bool_t   fExternalNmissed;                   ///< Set to true if want to give Nmissing from the MassStructureTask
  TRandom3 *fRndm;                             ///< TRandom3 object
  Int_t    fNMissedTracks;                     ///< Track missed in reconstruction calculated from external input (to be improved)
  //  -- now done in the analysis task from the difference of the integral of the pt distr of the constituents of the particle level jet and the reco level jet
  TClonesArray *fArrayTrackCorr;               ///< TClonesArray containing the jet constituents after correction
  TString  fPartArrayN;                        ///< Array of particles used for jet reconstruction at particle level (need to make it transient probably)
  //book-keeping object filled inside Eval()
  TProfile *fpAppliedEfficiency;               ///< Control profile efficiency
  THnSparse*fhNmissing;                        ///< pTjet vs number of added constituents (depends on settings) versus Nmissed constituents
  TList    *fListOfOutput;                     ///< list containing all histograms
 private:
  ClassDef(AliEmcalJetByJetCorrection, 8) // jet-by-jet correction class
};
#endif
