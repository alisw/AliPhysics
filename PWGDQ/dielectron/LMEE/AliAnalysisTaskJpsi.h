#ifndef AliAnalysisTaskJpsi_H
#define AliAnalysisTaskJpsi_H

//########################################################################
//#                                                                      #
//#   Class to produce charm cocktail (analysis task)                    #
//#                                                                      #
//#  Authors:                                                            #
//#   Raphaelle Bailhache, Uni Frankfurt / Raphaelle.Bailhache@cern.ch   #
//#   Carsten Klein, Uni Frankfurt / Carsten.Klein@cern.ch               #
//#   Jerome Jung, Uni Frankfurt / s8286523@uni-frankfurt.de             #
//#   Sebastian Scheid, Uni Frankfurt / s.scheid@cern.ch                 #
//#                                                                      #
//########################################################################

// local files
#include "AliCocktailSmearing.h"
#include "TF1.h"
class TTree;
class TParticle;
class AliAODMCParticle;
class TH1D;
class TH2D;
class TList;


class AliAnalysisTaskJpsi : public AliCocktailSmearing {
  
public:
  AliAnalysisTaskJpsi(); ///< default constructor probably needed for AnalysisManager or such...
  AliAnalysisTaskJpsi(const Char_t* name); ///< named constructor which also creates input and output objects.
  virtual ~AliAnalysisTaskJpsi();
  
  void Init();
  void ConnectTree(TTree *tree);
  void DetermineWeights();
  void Fill();

  TList *GetOutputList() {return fOutputList;};
  TList *GetOutputListLow() {return fOutputList_Low;};
  TList *GetOutputListHigh() {return fOutputList_High;};

  void SetOldTree(Bool_t oldTree)  {fOldTree = oldTree;}; 

  void SetYmin(Double_t Ymin) {fYmin = Ymin;};
  void SetYmax(Double_t Ymax) {fYmax = Ymax;};

  void SetPtRange(Double_t minpt, Double_t maxpt) {fPtCutRange[0] = minpt; fPtCutRange[1] = maxpt;};
  void SetEtaRange(Double_t mineta, Double_t maxeta) {fEtaCutRange[0] = mineta; fEtaCutRange[1] = maxeta;};

  void SetScaleFunction(TF1 *func) {fScaleFunction = func;};
  void SetScaleFunctionLow(TF1 *func) {fScaleFunctionLow = func;};
  void SetScaleFunctionHigh(TF1 *func) {fScaleFunctionHigh = func;};
  
private:
  void CreateHistos();
    
  
protected:
  TTree*              fTree;    //! tree
  TParticle*          fMotherOld;     //! mother
  TParticle*          fDaughter1Old;  //! daughter 1
  TParticle*          fDaughter2Old;  //! daughter 2
  AliAODMCParticle*       fMother;     //! mother
  AliAODMCParticle*       fDaughter1;  //! daughter 1
  AliAODMCParticle*       fDaughter2;  //! daughter 2

  Bool_t              fOldTree;    // Old format of the tree with TParticle
  
  Double_t            fYmin;       // Ymin for the normalization
  Double_t            fYmax;       // Ymax for the normalization

  Double_t            fPtCutRange[2];   // pt cut for electrons (Min and Max)
  Double_t            fEtaCutRange[2];  // eta cut for electrons (Min and Max)

  TF1*                fScaleFunction;  // scale function 1/N dN/dptdy
  TF1*                fScaleFunctionLow; // scale function 1/N dN/dptdy
  TF1*                fScaleFunctionHigh; // scale function 1/N dN/dptdy

  
  TH1D *fJPsiPt;                     // Histo for normalization
  TH2D *fMee_Ptee_Jpsi_rec_die;      // ptee,mee histos
  TH2D *fMee_Ptee_Jpsi_gen_die;
  TH2D *fMee_Ptee_Jpsi_rec_die_low;
  TH2D *fMee_Ptee_Jpsi_gen_die_low;
  TH2D *fMee_Ptee_Jpsi_rec_die_high;
  TH2D *fMee_Ptee_Jpsi_gen_die_high;
  TH2D *fMee_Ptee_Jpsi_rec_rad;
  TH2D *fMee_Ptee_Jpsi_gen_rad;
  TH2D *fMee_Ptee_Jpsi_rec_rad_low;
  TH2D *fMee_Ptee_Jpsi_gen_rad_low;
  TH2D *fMee_Ptee_Jpsi_rec_rad_high;
  TH2D *fMee_Ptee_Jpsi_gen_rad_high;

  TList       *fOutputList; //! Output list
  TList       *fOutputList_Low; //! Output list
  TList       *fOutputList_High; //! Output list

  AliAnalysisTaskJpsi(const AliAnalysisTaskJpsi &c); // not implemented
  AliAnalysisTaskJpsi& operator= (const AliAnalysisTaskJpsi &c); // not implemented
  
  ClassDef(AliAnalysisTaskJpsi,2)
};

#endif

