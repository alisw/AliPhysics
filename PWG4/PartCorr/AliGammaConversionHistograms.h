#ifndef ALIGAMMACONVERSIONHISTOGRAMS_H
#define ALIGAMMACONVERSIONHISTOGRAMS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

#include "TH1F.h"
#include "TH2F.h"
#include <Riostream.h>
#include <vector>
#include "TString.h"
#include "TList.h"

class AliGammaConversionHistograms{

 public: 
  
  AliGammaConversionHistograms();                                                         //constructor
  AliGammaConversionHistograms(const AliGammaConversionHistograms & g);                   //copy constructor
  AliGammaConversionHistograms & operator = (const AliGammaConversionHistograms & g);     //assignment operator
  virtual ~AliGammaConversionHistograms();                                                //virtual destructor
  

  TList * GetOutputContainer();
  
  Int_t GetRBin(Double_t radius);
  Int_t GetPhiBin(Double_t phi);

  //Setters/Initializers

  void Initialize_MappingValues(Int_t nPhiHistograms, Int_t nRHistograms, Int_t nBinsR, Double_t minRadius, Double_t maxRadius,Int_t nBinsPhi, Double_t minPhi, Double_t maxPhi);

  void Initialize_MC_EP_R(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_EP_Z_R(Int_t nXBins,Double_t firstX,Double_t lastX,Int_t nYBins,Double_t firstY,Double_t lastY,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_EP_X_Y(Int_t nXBins,Double_t firstX,Double_t lastX,Int_t nYBins,Double_t firstY,Double_t lastY,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_EP_OpeningAngle(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");

  void Initialize_MC_E_Energy(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_E_Pt(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_E_Eta(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_E_Phi(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");

  void Initialize_MC_P_Energy(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_P_Pt(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_P_Eta(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_P_Phi(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");

  void Initialize_MC_Gamma_Energy(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Gamma_Pt(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Gamma_Eta(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Gamma_Phi(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");

  void Initialize_MC_DirectGamma_Energy(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_DirectGamma_Pt(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_DirectGamma_Eta(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_DirectGamma_Phi(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");

  //mapping
  void Initialize_MappingHistograms(Int_t nPhiHistograms, Int_t nRHistograms,Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle);

  void Initialize_MC_Match_Gamma_Eta(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Match_Gamma_Phi(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Match_Gamma_Pt(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Match_Gamma_Energy(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Match_Gamma_Mass(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Match_Gamma_OpeningAngle(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Match_Gamma_R(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Match_Gamma_Z_R(Int_t nXBins,Double_t firstX,Double_t lastX,Int_t nYBins,Double_t firstY,Double_t lastY,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Match_Gamma_X_Y(Int_t nXBins,Double_t firstX,Double_t lastX,Int_t nYBins,Double_t firstY,Double_t lastY,TString xAxisTitle="", TString yAxisTitle="");

  void Initialize_MC_Pi0_Eta(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Pi0_Phi(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Pi0_Pt(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Pi0_Energy(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Pi0_Mass(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Pi0_OpeningAngleGamma(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Pi0_R(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Pi0_Z_R(Int_t nXBins,Double_t firstX,Double_t lastX,Int_t nYBins,Double_t firstY,Double_t lastY,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Pi0_X_Y(Int_t nXBins,Double_t firstX,Double_t lastX,Int_t nYBins,Double_t firstY,Double_t lastY,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Pi0Secondaries_X_Y(Int_t nXBins,Double_t firstX,Double_t lastX,Int_t nYBins,Double_t firstY,Double_t lastY,TString xAxisTitle="", TString yAxisTitle="");

  void Initialize_MC_Eta_Eta(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Eta_Phi(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Eta_Pt(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Eta_Energy(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Eta_Mass(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Eta_OpeningAngleGamma(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Eta_R(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Eta_Z_R(Int_t nXBins,Double_t firstX,Double_t lastX,Int_t nYBins,Double_t firstY,Double_t lastY,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_MC_Eta_X_Y(Int_t nXBins,Double_t firstX,Double_t lastX,Int_t nYBins,Double_t firstY,Double_t lastY,TString xAxisTitle="", TString yAxisTitle="");
    
  // esd

  void Initialize_ESD_EP_R(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_EP_Z_R(Int_t nXBins,Double_t firstX,Double_t lastX,Int_t nYBins,Double_t firstY,Double_t lastY,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_EP_X_Y(Int_t nXBins,Double_t firstX,Double_t lastX,Int_t nYBins,Double_t firstY,Double_t lastY,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_EP_OpeningAngle(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");

  void Initialize_ESD_E_Energy(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_E_Pt(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_E_Eta(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_E_Phi(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");

  void Initialize_ESD_P_Energy(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_P_Pt(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_P_Eta(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_P_Phi(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");


  void Initialize_ESD_Gamma_Energy(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Gamma_Pt(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Gamma_Eta(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Gamma_Phi(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");

  void Initialize_ESD_Match_Gamma_OpeningAngle(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Match_Gamma_Energy(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Match_Gamma_Pt(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Match_Gamma_Eta(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Match_Gamma_Phi(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Match_Gamma_Mass(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Match_Gamma_Width(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Match_Gamma_Chi2(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Match_Gamma_NDF(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Match_Gamma_R(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Match_Gamma_Z_R(Int_t nXBins,Double_t firstX,Double_t lastX,Int_t nYBins,Double_t firstY,Double_t lastY,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Match_Gamma_X_Y(Int_t nXBins,Double_t firstX,Double_t lastX,Int_t nYBins,Double_t firstY,Double_t lastY,TString xAxisTitle="", TString yAxisTitle="");


  void Initialize_ESD_Pi0_OpeningAngleGamma(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Pi0_Energy(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Pi0_Pt(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Pi0_Eta(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Pi0_Phi(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Pi0_Mass(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Pi0_R(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Pi0_Z_R(Int_t nXBins,Double_t firstX,Double_t lastX,Int_t nYBins,Double_t firstY,Double_t lastY,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Pi0_X_Y(Int_t nXBins,Double_t firstX,Double_t lastX,Int_t nYBins,Double_t firstY,Double_t lastY,TString xAxisTitle="", TString yAxisTitle="");


  void Initialize_ESD_Eta_OpeningAngleGamma(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Eta_Energy(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Eta_Pt(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Eta_Eta(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Eta_Phi(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Eta_Mass(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Eta_R(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Eta_Z_R(Int_t nXBins,Double_t firstX,Double_t lastX,Int_t nYBins,Double_t firstY,Double_t lastY,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Eta_X_Y(Int_t nXBins,Double_t firstX,Double_t lastX,Int_t nYBins,Double_t firstY,Double_t lastY,TString xAxisTitle="", TString yAxisTitle="");

  void Initialize_ESD_Background_OpeningAngleGamma(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Background_Energy(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Background_Pt(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Background_Eta(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Background_Phi(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Background_Mass(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Background_R(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Background_Z_R(Int_t nXBins,Double_t firstX,Double_t lastX,Int_t nYBins,Double_t firstY,Double_t lastY,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_ESD_Background_X_Y(Int_t nXBins,Double_t firstX,Double_t lastX,Int_t nYBins,Double_t firstY,Double_t lastY,TString xAxisTitle="", TString yAxisTitle="");


  void Initialize_Resolution_dPt(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle);
  void Initialize_Resolution_dR(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle);
  void Initialize_Resolution_dZ(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle);
  void Initialize_Resolution_dR_dPt(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle);

  void Initialize_Resolution_MC_Pt(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_Resolution_MC_R(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_Resolution_MC_Z(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");

  void Initialize_Resolution_ESD_Pt(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_Resolution_ESD_R(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_Resolution_ESD_Z(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");

  void Initialize_NumberOfV0s(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_NumberOfSurvivingV0s(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");

  void Initialize_V0MassDebugCut1(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_V0MassDebugCut2(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_V0MassDebugCut3(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_V0MassDebugCut4(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_V0MassDebugCut5(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_V0MassDebugCut6(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_V0MassDebugCut7(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");
  void Initialize_V0MassDebugCut8(Int_t nXBins,Double_t firstX,Double_t lastX,TString xAxisTitle="", TString yAxisTitle="");

 private:
  TList* fOutputContainer;

  Int_t fNPhiIndex;
  Int_t fNRIndex;
  Double_t fMinRadius;
  Double_t fMaxRadius;
  Double_t fDeltaR;
  Double_t fMinPhi;
  Double_t fMaxPhi;
  Double_t fDeltaPhi;
  

  
  // Pure MonteCarlo histograms
  TH1F * fMC_EP_R;                    //! transient
  TH2F * fMC_EP_Z_R;                  //! transient 
  TH2F * fMC_EP_X_Y;                  //! transient
  TH1F * fMC_EP_OpeningAngle;         //! transient

  TH1F * fMC_E_Energy;                //! transient
  TH1F * fMC_E_Pt;                    //! transient
  TH1F * fMC_E_Eta;                   //! transient
  TH1F * fMC_E_Phi;                   //! transient

  TH1F * fMC_P_Energy;                //! transient
  TH1F * fMC_P_Pt;                    //! transient
  TH1F * fMC_P_Eta;                   //! transient
  TH1F * fMC_P_Phi;                   //! transient

  TH1F * fMC_Gamma_Energy;            //! transient
  TH1F * fMC_Gamma_Pt;                //! transient
  TH1F * fMC_Gamma_Eta;               //! transient
  TH1F * fMC_Gamma_Phi;               //! transient

  TH1F * fMC_DirectGamma_Energy;      //! transient
  TH1F * fMC_DirectGamma_Pt;          //! transient
  TH1F * fMC_DirectGamma_Eta;         //! transient
  TH1F * fMC_DirectGamma_Phi;         //! transient

  //mapping
  //begin double vector
  typedef vector<TH2F *> AliConversionMappingVector; //! transient  
  vector<AliConversionMappingVector> fMC_Mapping;    //! transient
  //end double vector

  vector<TH2F *> fMC_Mapping_Phi;     //! transient
  vector<TH2F *> fMC_Mapping_R;       //! transient

  TH1F * fMC_Match_Gamma_Eta;         //! transient
  TH1F * fMC_Match_Gamma_Phi;         //! transient
  TH1F * fMC_Match_Gamma_Pt;          //! transient
  TH1F * fMC_Match_Gamma_Energy;      //! transient
  TH1F * fMC_Match_Gamma_Mass;        //! transient
  TH1F * fMC_Match_Gamma_OpeningAngle;//! transient
  TH1F * fMC_Match_Gamma_R;           //! transient
  TH2F * fMC_Match_Gamma_Z_R;         //! transient
  TH2F * fMC_Match_Gamma_X_Y;         //! transient


  TH1F * fMC_Pi0_Eta;                 //! transient
  TH1F * fMC_Pi0_Phi;                 //! transient
  TH1F * fMC_Pi0_Pt;                  //! transient
  TH1F * fMC_Pi0_Energy;              //! transient
  TH1F * fMC_Pi0_Mass;                //! transient Not filled, no point, we know the montecarlo mass
  TH1F * fMC_Pi0_OpeningAngleGamma;   //! transient
  TH1F * fMC_Pi0_R;                   //! transient
  TH2F * fMC_Pi0_Z_R;                 //! transient
  TH2F * fMC_Pi0_X_Y;                 //! transient

  TH1F * fMC_Pi0Secondaries_Eta;                 //! transient
  TH1F * fMC_Pi0Secondaries_Phi;                 //! transient
  TH1F * fMC_Pi0Secondaries_Pt;                  //! transient
  TH1F * fMC_Pi0Secondaries_Energy;              //! transient
  TH1F * fMC_Pi0Secondaries_Mass;                //! transient Not filled, no point, we know the monteacrlo mass
  TH1F * fMC_Pi0Secondaries_OpeningAngleGamma;   //! transient
  TH1F * fMC_Pi0Secondaries_R;                   //! transient
  TH2F * fMC_Pi0Secondaries_Z_R;                 //! transient
  TH2F * fMC_Pi0Secondaries_X_Y;

  TH1F * fMC_Eta_Eta;                            //! transient
  TH1F * fMC_Eta_Phi;                            //! transient
  TH1F * fMC_Eta_Pt;                             //! transient
  TH1F * fMC_Eta_Energy;                         //! transient
  TH1F * fMC_Eta_Mass;                           //! transient Not Filled, no point we know the montecarlo mass
  TH1F * fMC_Eta_OpeningAngleGamma;              //! transient 
  TH1F * fMC_Eta_R;                              //! transient We have very few eta secondaries, so the question is if we keep this
  TH2F * fMC_Eta_Z_R;                            //! transient Same here, do we really need it?
  TH2F * fMC_Eta_X_Y;                            //! transient all the etas has their vertex in a square in the collision point
    
  // Histograms from esd tracks
  TH1F * fESD_EP_R;                              //! transient  
  TH2F * fESD_EP_Z_R;                            //! transient
  TH2F * fESD_EP_X_Y;                            //! transient
  TH1F * fESD_EP_OpeningAngle;                   //! transient

  TH1F * fESD_E_Energy;                          //! transient
  TH1F * fESD_E_Pt;                              //! transient
  TH1F * fESD_E_Eta;                             //! transient
  TH1F * fESD_E_Phi;                             //! transient

  TH1F * fESD_P_Energy;                          //! transient
  TH1F * fESD_P_Pt;                              //! transient
  TH1F * fESD_P_Eta;                             //! transient
  TH1F * fESD_P_Phi;                             //! transient


  TH1F * fESD_Gamma_Energy;                      //! transient
  TH1F * fESD_Gamma_Pt;                          //! transient
  TH1F * fESD_Gamma_Eta;                         //! transient
  TH1F * fESD_Gamma_Phi;                         //! transient


  //mapping
  //begin double vector
  typedef vector<TH2F *> fESDPhiRVector; //! transient  
  vector<fESDPhiRVector> fESD_Mapping;    //! transient
  //end double vector
  vector<TH2F *> fESD_Mapping_Phi;               //! transient
  vector<TH2F *> fESD_Mapping_R;                 //! transient

  TH1F * fESD_Match_Gamma_OpeningAngle;          //! transient 
  TH1F * fESD_Match_Gamma_Energy;                //! transient
  TH1F * fESD_Match_Gamma_Pt;                    //! transient
  TH1F * fESD_Match_Gamma_Eta;                   //! transient
  TH1F * fESD_Match_Gamma_Phi;                   //! transient
  TH1F * fESD_Match_Gamma_Mass;                  //! transient
  TH1F * fESD_Match_Gamma_Width;                 //! transient
  TH1F * fESD_Match_Gamma_Chi2;                  //! transient
  TH1F * fESD_Match_Gamma_NDF;                   //! transient
  TH1F * fESD_Match_Gamma_R;                     //! transient
  TH2F * fESD_Match_Gamma_Z_R;                   //! transient
  TH2F * fESD_Match_Gamma_X_Y;                   //! transient


  TH1F * fESD_Pi0_OpeningAngleGamma;             //! transient
  TH1F * fESD_Pi0_Energy;                        //! transient
  TH1F * fESD_Pi0_Pt;                            //! transient
  TH1F * fESD_Pi0_Eta;                           //! transient
  TH1F * fESD_Pi0_Phi;                           //! transient
  TH1F * fESD_Pi0_Mass;                          //! transient
  TH1F * fESD_Pi0_R;                             //! transient
  TH2F * fESD_Pi0_Z_R;                           //! transient
  TH2F * fESD_Pi0_X_Y;                           //! transient

  TH1F * fESD_Eta_OpeningAngleGamma;             //! transient
  TH1F * fESD_Eta_Energy;                        //! transient
  TH1F * fESD_Eta_Pt;                            //! transient
  TH1F * fESD_Eta_Eta;                           //! transient
  TH1F * fESD_Eta_Phi;                           //! transient
  TH1F * fESD_Eta_Mass;                          //! transient
  TH1F * fESD_Eta_R;                             //! transient
  TH2F * fESD_Eta_Z_R;                           //! transient
  TH2F * fESD_Eta_X_Y;                           //! transient

  TH1F * fESD_Background_OpeningAngleGamma;      //! transient
  TH1F * fESD_Background_Energy;                 //! transient
  TH1F * fESD_Background_Pt;                     //! transient
  TH1F * fESD_Background_Eta;                    //! transient
  TH1F * fESD_Background_Phi;                    //! transient
  TH1F * fESD_Background_Mass;                   //! transient
  TH1F * fESD_Background_R;                      //! transient
  TH2F * fESD_Background_Z_R;                    //! transient
  TH2F * fESD_Background_X_Y;                    //! transient

  TH2F * fResolution_dPt;                        //! transient
  TH2F * fResolution_dR;                         //! transient
  TH2F * fResolution_dZ;                         //! transient
  
  TH2F * fResolution_dR_dPt;                     //! transient

  TH1F * fResolution_MC_Pt;                      //! transient
  TH1F * fResolution_MC_R;                       //! transient
  TH1F * fResolution_MC_Z;                       //! transient

  TH1F * fResolution_ESD_Pt;                     //! transient
  TH1F * fResolution_ESD_R;                      //! transient
  TH1F * fResolution_ESD_Z;                      //! transient

  TH1F * fNumberOfV0s;                           //! transient
  TH1F * fNumberOfSurvivingV0s;                  //! transient

  //  debug histograms
  TH1F * fV0MassDebugCut1;                       //! transient
  TH1F * fV0MassDebugCut2;                       //! transient
  TH1F * fV0MassDebugCut3;                       //! transient
  TH1F * fV0MassDebugCut4;                       //! transient
  TH1F * fV0MassDebugCut5;                       //! transient
  TH1F * fV0MassDebugCut6;                       //! transient
  TH1F * fV0MassDebugCut7;                       //! transient
  TH1F * fV0MassDebugCut8;                       //! transient


  friend class AliAnalysisTaskGammaConversion;
  friend class AliV0Reader;

  ClassDef(AliGammaConversionHistograms,0)
} ;


#endif



