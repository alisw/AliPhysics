#ifndef ALIANALYSISKINKESDMC_H
#define ALIANALYSISKINKESDMC_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisKinkESDMC class
//         This task is an example of an analysis task
//                  for kink topology Study
//          Authors: Martha Spyropoulou-Stassinaki
//           and members of the Greek group at the
//          Physics Department of Athens University
//                    mspyrop@phys.uoa.gr
//-----------------------------------------------------------------

class AliESDVertex;
class AliESDEvent;
class TTree;
class TH1F;
class TH2F;
class TH1D;
class TH2D;
class AliPID;
class AliAnalysisTaskSE;

class AliAnalysisKinkESDMC : public AliAnalysisTaskSE {
 public:
  AliAnalysisKinkESDMC();
  AliAnalysisKinkESDMC(const char *name);
  virtual ~AliAnalysisKinkESDMC() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
 
  Float_t  GetSigmaToVertex(AliESDtrack* esdTrack) const;
  const AliESDVertex *GetEventVertex(const AliESDEvent* esd) const;
  
 private:
   AliESDEvent *fESD;    //!ESD object
   TH1F        *fHistPtESD; //!Pt spectrum of all ESD inside eta, Pt cuts
   TH1F        *fHistPt; //!Pt spectrum of all ESD tracks
   TH1F        *fHistQtAll; //!Qt spectrum of all kinks
   TH1F        *fHistQt1; //!Qt spectrum of Kaon selected sample
   TH1F        *fHistQt2; //!Qt spectrum in Qt region of kaons
   TH1F        *fHistPtKaon; //!Pt Kaon spectrum of clean sample
   TH1F        *fHistPtKPDG; //!Pt Kaon spectrum , confirmed by  PDG,inside kaon Qt region
   TH1F        *fHistEta; //!Eta spectrum of all kinks
   TH1F        *fHistEtaK; //!Eta spectrum of kaons selected by kink topology
   TH1F        *fptKMC; //!Pt Kaon spectrum MC, inside eta and pt cuts 
   TH1F        *fMultiplMC; //!charge multipl MC 
   TH1F        *fESDMult; //!ESD charged mult
   TH1F        *fgenpt; //!Pt Kaon-Kink->mu  spectrum , MC, inside eta, Pt, radius cuts
   TH1F        *frad; //!radius of kinks,  MC , inside the eta nad Pt cuts 
   TH1F        *fKinkKaon; //!Pt of PDG Kaons inside the selcted ones by the KInk topology 
   TH1F        *fKinkKaonBg; //!Pt of the BG inside the kink-Kaon identified spectrum
   TH1F        *fM1kaon; //!inv mass of kink-tracks taken as kaons decaying to  mu + neutrino
   TH1F        *fgenPtEtR; //!MC Pt spectrum of kaons decaying to muon+neutrino and pi +pi, inside eta,Pt,Rad cuts
   TH1F        *fPtKink; //!Pt  spectrum   of all kinks  from track bank
   TH1F        *fptKink; //!Pt  spectrum of all kinks from kink bank
   TH2F        *fcodeH ; //!PDG code(mother)  vrs PDG dcode(daughter) of kinks with Qt <0.12 (fake)
   TH2F        *fdcodeH ; //!Kinks, code  vrs dcode of BG,if mother code is 321 and daughter code > 
   TH2F        *fAngMomK; //! Decay angle vrs Mother Mom for pdg kaons
   TH2F        *fAngMomPi; //! Decay angle vrs Mother Mom for pdg pions
   TH1D        *fRpr; //! Radius of VTX at Y , X plane              
   TH1D        *fZpr; //! Z distrio of main vertex                  
   TH2F        *fAngMomKC; //! Decay angle vrs Mother Mom for pdg kaons, inside the selected sample
   TH2F        *fZvXv; //! two dime of Z vrs X of vtx main           
   TH2F        *fZvYv; //! two dime of Z vrs Y of vtx main           
   TH2F        *fXvYv; //! two dime of X vrs Y of main tracks vtx main           
   TH1F        *fPtPrKink; //! pt of Primary PDG kaons inside the selected ones by the kink topology              
   TF1         *f1;
   TF1         *f2;
  TList        *fListOfHistos; //! list of histos

  AliAnalysisKinkESDMC(const AliAnalysisKinkESDMC&); // not implemented
  AliAnalysisKinkESDMC& operator=(const AliAnalysisKinkESDMC&); // not implemented

  ClassDef(AliAnalysisKinkESDMC, 1); // example of analysis
};

#endif
