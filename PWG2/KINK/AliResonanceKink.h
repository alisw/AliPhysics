#ifndef ALIRESONANCEKINK_H
#define ALIRESONANCEKINK_H

/*  See cxx source for full Copyright notice */

//------------------------------------------------------------------------------
//                   class AliResonanceKink
//         This task is an example of an analysis task
//        for analysing resonances having one kaon kink
//Author: Paraskevi Ganoti, University of Athens (pganoti@phys.uoa.gr)
//------------------------------------------------------------------------------

class TF1;
class AliESDEvent;
class AliESDtrack;
class AliESDVertex;
class AliMCEvent;
class TList;
class TString;

class AliResonanceKink : public TObject {
 public:
 
  enum ResonanceType {kPhi=333, kKstar0=313, kLambda1520=3124};
  enum DaughterType {kdaughterPion=211, kdaughterKaon=321, kdaughterProton=2212};
  
  AliResonanceKink();
  AliResonanceKink(Int_t nbins, Float_t nlowx, Float_t nhighx, Int_t daughter1, Int_t daughter2, Int_t resonancePDG);
  virtual ~AliResonanceKink();
  
  TList* GetHistogramList();
  void Analyse(AliESDEvent* esd, AliMCEvent* mcEvent);  
  Float_t GetSigmaToVertex(AliESDtrack* esdTrack) const ; 
  const AliESDVertex *GetEventVertex(const AliESDEvent* esd) const;
  void SetAnalysisType(TString type) {fAnalysisType=type;}
  void SetPDGCodes(Int_t d1, Int_t d2, Int_t res) {fdaughter1pdg=d1; fdaughter2pdg=d2; fresonancePDGcode=res;}
  void InitOutputHistograms(Int_t nbins, Float_t nlowx, Float_t nhighx);
 
 private:
 
  TList       *fListOfHistos; //! List 
  TH1D        *fOpeningAngle; //! Opening  
  TH1D        *fInvariantMass; //! invMass spectrum   
  TH1D        *fInvMassTrue; //! invMassTrue spectrum  
  TH1D        *fPhiBothKinks; //! bothKaonsKinks   
  TH1D        *fRecPt; //! pT spectrum  
  TH1D        *fRecEta; //! Eta spectrum
  TH2D        *fRecEtaPt; //! Eta pT spectrum  
  TH1D        *fSimPt; //! pT Sim spectrum  
  TH1D        *fSimEta; //! Eta Sim spectrum
  TH2D        *fSimEtaPt; //! Eta pT Sim spectrum 
  TH1D        *fSimPtKink; //! pT Sim one kaon kink spectrum  
  TH1D        *fSimEtaKink; //! Eta Sim one kaon kink spectrum spectrum
  TH2D        *fSimEtaPtKink; //! Eta pT Sim one kaon kink spectrum   
  TH1D        *fhdr ; //! radial impact  
  TH1D        *fhdz ; //! z impact
  TF1         *f1;
  TF1         *f2;  
  TString     fAnalysisType;//"ESD" or "MC"
  TH1D        *fvtxz ; //! vtx z component
  Int_t       fNbins; //! bins
  Float_t     fLowX;//! lowx
  Float_t     fHighX; //! high x
  Int_t       fdaughter1pdg;
  Int_t       fdaughter2pdg;  
  Int_t       fresonancePDGcode;
  
  AliResonanceKink(const AliResonanceKink&); // not implemented
  AliResonanceKink& operator=(const AliResonanceKink&); // not implemented

  ClassDef(AliResonanceKink, 1); // example of analysis
};

#endif
