#ifndef AliAnalysisTaskEmcalJetSpectra_h
#define AliAnalysisTaskEmcalJetSpectra_h

// $Id$


class TH1F;
class TH2F;
class THnSparse;
class AliLocalRhoParameter;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskEmcalJetSpectra : public AliAnalysisTaskEmcalJet {
 public:
  AliAnalysisTaskEmcalJetSpectra();
  AliAnalysisTaskEmcalJetSpectra(const char *name);
  virtual ~AliAnalysisTaskEmcalJetSpectra() {}
 
  
  virtual void           UserCreateOutputObjects();

 protected:
  Bool_t                 Run();
  virtual Int_t          GetCentBin(Double_t cent) const;
  Float_t                RelativePhi(Double_t mphi,Double_t vphi) const;
  Float_t                RelativeEPJET(Double_t jetAng, Double_t EPAng) const;
  Double_t	         fLocalRhoVal;

 private:
  TH2F                  *fHistRhovsCent; //!
  TH2F                  *fHistNjetvsCent;          //!number of jets versus Centrality
  TH2F                  *fHistJetPtvsTrackPt[6];//!
  TH2F                  *fHistRawJetPtvsTrackPt[6];//!
  TH1F                  *fHistTrackPt[6];//!
  TH1F                  *fHistEP0[6];//!
  TH1F                  *fHistEP0A[6];//!
  TH1F                  *fHistEP0C[6];//!
  TH2F                  *fHistEPAvsC[6];//!
  TH2F                  *fHistJetPtvsdEP[6];//!
  TH2F                  *fHistJetPtvsdEPBias[6];//!
  TH2F                  *fHistJetPtvsEP[6];//!
  TH2F                  *fHistJetPtvsEPBias[6];//!
  TH2F                  *fHistRhovsEP[6]; //!
  TH1F                  *fHistCorJetPtfromLocalRho[6]; //!
  TH1F                  *fHistCorJetPtfromGlobalRho[6]; //!

  TH2F                  *fHistGLvsLOCrho; //!         // Global vs Local Rho distribution
  TH2F                  *fHistRhovsdEPLOC; //! 
  TH2F                  *fHistRhovsdEPGL; //! 
  TH2F                  *fHistJetPtvsdEPLOC; //! 
  TH2F                  *fHistJetPtvsdEPGL; //! 
  TH2F                  *fHistRhovsEPLOC; //! 
  TH2F                  *fHistRhovsEPGL; //! 
  TH2F                  *fHistJetPtvsEPLOC; //!  
  TH2F                  *fHistJetPtvsEPGL; //! 
  TH1F                  *fHistCorJetPt;  //!            // (Njets) vs Corrected Jet Pt (local rho)
  TH1F                  *fHistCorJetPtGL; //!           // (Njets) vs Corrected Jet Pt (global rho)

  TH1F                  *fHistCorJetPtfromLocalRhoIN[6]; //! 
  TH1F                  *fHistCorJetPtfromLocalRhoOUT[6]; //! 
  TH1F                  *fHistCorJetPtfromGlobalRhoIN[6]; //! 
  TH1F                  *fHistCorJetPtfromGlobalRhoOUT[6]; //! 

  TH2F                  *fHistRhodEPcentLOC[6]; //! 
  TH2F                  *fHistRhodEPcentGL[6]; //! 
  TH2F                  *fHistCorJetPtdEPcentLOC[6]; //! 
  TH2F                  *fHistCorJetPtdEPcentGL[6]; //! 

  TH2F                  *fHistRhoEPcentLOC[6]; //! 
  TH2F                  *fHistRhoEPcentGL[6]; //! 
  TH2F                  *fHistCorJetPtEPcentLOC[6]; //! 
  TH2F                  *fHistCorJetPtEPcentGL[6]; //! 


  AliAnalysisTaskEmcalJetSpectra(const AliAnalysisTaskEmcalJetSpectra&); // not implemented
  AliAnalysisTaskEmcalJetSpectra& operator=(const AliAnalysisTaskEmcalJetSpectra&); // not implemented
  
  ClassDef(AliAnalysisTaskEmcalJetSpectra, 5); // Emcal jet spectra task
};
#endif
