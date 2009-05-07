#ifndef ALIRESONANCEKINKPID_H
#define ALIRESONANCEKINKPID_H

/*  See cxx source for full Copyright notice */

//------------------------------------------------------------------------------
//                   class AliResonanceKinkPID
//         This task is an example of an analysis task
//        for analysing resonances having one kaon kink
//Author: Paraskevi Ganoti, University of Athens (pganoti@phys.uoa.gr)
//------------------------------------------------------------------------------

class TF1;
class TString;
class TTree;
class AliESDEvent;
class AliESDtrack;

class AliResonanceKinkPID : public AliAnalysisTaskSE {
 public:
  AliResonanceKinkPID();
  AliResonanceKinkPID(const char *name);
  virtual ~AliResonanceKinkPID() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  Float_t GetSigmaToVertex(AliESDtrack* esdTrack) const ; 
  const AliESDVertex *GetEventVertex(const AliESDEvent* esd) const;
  
 private:
  AliESDEvent *fESD;    //! ESD object
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
  
  AliResonanceKinkPID(const AliResonanceKinkPID&); // not implemented
  AliResonanceKinkPID& operator=(const AliResonanceKinkPID&); // not implemented

  ClassDef(AliResonanceKinkPID, 1); // example of analysis
};

#endif
