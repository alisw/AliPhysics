#ifndef AliAnalysisTaskEmcalJetHMEC_cxx
#define AliAnalysisTaskEmcalJetHMEC_cxx


class TList;
class TH1;
class TH2;
class THnSparse;
class AliESDEvent;
class AliEventPoolManager;

#include "AliLog.h"

#include "AliAnalysisTaskSE.h"
//#include "/project/projectdirs/alice/tschuste/AliceSoftware/aliroot/train_jet/ANALYSIS/AliAnalysisTaskSE.h"

class AliAnalysisTaskEmcalJetHMEC : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEmcalJetHMEC();
  AliAnalysisTaskEmcalJetHMEC(const char *name);
  virtual ~AliAnalysisTaskEmcalJetHMEC() {}
  
  virtual void      UserCreateOutputObjects();
  virtual Double_t   RelativePhi(Double_t mphi,Double_t vphi);
  virtual void      UserExec(Option_t *option);
  virtual void      Terminate(Option_t *);
   virtual THnSparse* NewTHnSparseF(const char* name, UInt_t entries);
   virtual void       GetDimParams(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);

  virtual void      SetTracksName(const char *n) {fTracksName=n;}
  virtual void      SetJetsName(const char *jn) {fJetsName=jn;}

  virtual void           SetAreaCut(Double_t a)                   { fAreacut    = a; }
  virtual void           SetJetEta(Double_t emin, Double_t emax)  { fEtamin = emin; fEtamax = emax; }
  virtual void           SetJetPhi(Double_t pmin, Double_t pmax)  { fPhimin = pmin; fPhimax = pmax; }
   virtual void     SetEventMixing(Int_t yesno){fDoEventMixing=yesno;}
    virtual	void    SetMixingTracks(Int_t tracks) { fMixingTracks = tracks; }





 protected:
  virtual Int_t          GetCentBin(Double_t cent) const;
  virtual Int_t          GetEtaBin(Double_t eta) const;
  virtual Int_t          GetpTjetBin(Double_t pt) const;

 private:
  TString      fTracksName;  //name of tracks collection
  TString      fJetsName;  //name of Jet collection
  Double_t               fPhimin;                  // phi min
  Double_t               fPhimax;                  // phi max
  Double_t               fEtamin;                  // eta min
  Double_t               fEtamax;                  // eta max
  Double_t               fAreacut;                 // area cut
   Int_t   fDoEventMixing;
    Int_t  		fMixingTracks;		// size of track buffer for event mixing
    TObjArray* CloneAndReduceTrackList(TObjArray* tracks);


  AliESDEvent *fESD;    //! ESD object
  AliEventPoolManager   *fPoolMgr;
  TList       *fOutputList; //! Output list
  TH1        *fHistTrackPt; //! Pt spectrum
  TH1         *fHistCentrality;
  TH2         *fHistJetEtaPhi;
  TH2         *fHistTrackEtaPhi;
  TH2         *fHistJetHEtaPhi;
   Int_t   fNevents;          // number of events
   Int_t   fTindex;           // index reference
   Int_t   fTrigBufferIndex;  //index for the buffering
   Int_t   fCountAgain;       //index for the buffering

  TH1         *fHistJetPt[6];
  TH1         *fHistJetPtBias[6];
  TH1         *fHistJetPtTT[6];
  TH2         *fHistJetH[6][5][3];
  TH2         *fHistJetHBias[6][5][3];
  TH2         *fHistJetHTT[6][5][3];
   THnSparse *fhnMixedEvents;                //!mixed events matrix
   Double_t            fTrigBuffer[10][7];      //!buffer for triggers   

   
  AliAnalysisTaskEmcalJetHMEC(const AliAnalysisTaskEmcalJetHMEC&); // not implemented
  AliAnalysisTaskEmcalJetHMEC& operator=(const AliAnalysisTaskEmcalJetHMEC&); // not implemented
  
  ClassDef(AliAnalysisTaskEmcalJetHMEC, 4); 
};


class AliDPhiBasicParticleMEC : public AliVParticle
{
  public:
    AliDPhiBasicParticleMEC(Float_t eta, Float_t phi, Float_t pt, Short_t charge)
      : fEta(eta), fPhi(phi), fpT(pt), fCharge(charge)
    {
    }
    ~AliDPhiBasicParticleMEC() {}
    
    // kinematics
    virtual Double_t Px() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Py() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Pz() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Pt() const { return fpT; }
    virtual Double_t P() const { AliFatal("Not implemented"); return 0; }
    virtual Bool_t   PxPyPz(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

    virtual Double_t Xv() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Yv() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Zv() const { AliFatal("Not implemented"); return 0; }
    virtual Bool_t   XvYvZv(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

    virtual Double_t OneOverPt()  const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Phi()        const { return fPhi; }
    virtual Double_t Theta()      const { AliFatal("Not implemented"); return 0; }


    virtual Double_t E()          const { AliFatal("Not implemented"); return 0; }
    virtual Double_t M()          const { AliFatal("Not implemented"); return 0; }
    
    virtual Double_t Eta()        const { return fEta; }
    virtual Double_t Y()          const { AliFatal("Not implemented"); return 0; }
    
    virtual Short_t Charge()      const { return fCharge; }
    virtual Int_t   GetLabel()    const { AliFatal("Not implemented"); return 0; }
    // PID
    virtual Int_t   PdgCode()     const { AliFatal("Not implemented"); return 0; }      
    virtual const Double_t *PID() const { AliFatal("Not implemented"); return 0; }
    
  private:
    Float_t fEta;      // eta
    Float_t fPhi;      // phi
    Float_t fpT;       // pT
    Short_t fCharge;   // charge
    
    ClassDef( AliDPhiBasicParticleMEC, 1); // class which contains only quantities requires for this analysis to reduce memory consumption for event mixing
};


#endif
