#ifndef AliAnaVZEROEPFlatenning_cxx
#define AliAnaVZEROEPFlatenning_cxx

class AliVEvent;

#include "AliAnalysisTaskSE.h"

class AliAnaVZEROEPFlatenning : public AliAnalysisTaskSE {
 public:
  AliAnaVZEROEPFlatenning();
  AliAnaVZEROEPFlatenning(const char *name);
  virtual ~AliAnaVZEROEPFlatenning() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  virtual void Init();

  void SetMBTrigName(const char *name = "CPBI") {fMBTrigName = name;}
  void SetUsePhysSel(Bool_t usePhysSel) {fUsePhysSel = usePhysSel;}
  void SetInput(const char *filename);

  Double_t CalculateVZEROEventPlane(const AliVEvent *  event, Int_t ring, Float_t centrality, Double_t &qxTierce, Double_t &qyTierce) const;

 private:
  AliVEvent   *fEvent;    //! ESD ot AOD object
  TList       *fOutputList; //! Output list

  TString fMBTrigName; // MB trigger name (for evt sel)
  Bool_t  fUsePhysSel; // Use or not phys sel

  TProfile *fX2[11]; //! Profile histogram for Q^2_x
  TProfile *fY2[11]; //! Profile histogram for Q^2_y
  TProfile *fX2Y2[11]; //! Profile histogram for Q^2_x*Q^2_y
  TProfile *fCos8Psi[11]; //! Profile histogram for Cos(8*Psi)
  TProfile *fC2[8]; //! Profile histogram for Cos(2*phi)
  TProfile *fS2[8]; //! Profile histogram for Sin(2*phi)
  TProfile *fC4[8]; //! Profile histogram for Cos(4*phi)
  TProfile *fS4[8]; //! Profile histogram for Sin(4*phi)
  TProfile *fX2Corr[11]; //! Profile histogram for Q^2_x
  TProfile *fY2Corr[11]; //! Profile histogram for Q^2_y
  TProfile *fX2Y2Corr[11]; //! Profile histogram for Q^2_x*Q^2_y

  TProfile *fX2In[8]; // Profile histogram for Q^2_x (read from input file)
  TProfile *fY2In[8]; // Profile histogram for Q^2_y (read from input file)
  TProfile *fX2Y2In[8]; // Profile histogram for Q^2_x*Q^2_y (read from input file)
  TProfile *fCos8PsiIn[8]; // Profile histogram for Cos(8*Psi) (read from input file)

  TH2F *fPsiRingRawCentr[8]; //! Raw VZERO event plane vs centrality (ring-by-ring)
  TH2F *fPsiRingFlatCentr[8]; //! Flatenned with corrections on cumulants VZERO event plane vs centrality (ring-by-ring)
  TH2F *fPsiRingFlatFinalCentr[8]; //! Flatenned with corrections on cumulants and fourier VZERO event plane vs centrality (ring-by-ring)
  TH2F *fPsiARawCentr; //! Raw VZEROA event plane vs centrality
  TH2F *fPsiAFlatCentr; //! Flatenned with corrections on cumulants VZEROA event plane vs centrality 
  TH2F *fPsiCRawCentr; //! Raw VZEROC event plane vs centrality
  TH2F *fPsiCFlatCentr; //! Flatenned with corrections on cumulants VZEROC event plane vs centrality 
  TH2F *fPsiACRawCentr; //! Raw VZERO event plane vs centrality
  TH2F *fPsiACFlatCentr; //! Flatenned with corrections on cumulants VZERO event plane vs centrality 

  AliAnaVZEROEPFlatenning(const AliAnaVZEROEPFlatenning&); // not implemented
  AliAnaVZEROEPFlatenning& operator=(const AliAnaVZEROEPFlatenning&); // not implemented
  
  ClassDef(AliAnaVZEROEPFlatenning, 1); // VZERO analysis task for extraction EP flatenning params
};

#endif
