#ifndef ALIITSSIMULATIONSPDDUBNA_H
#define ALIITSSIMULATIONSPDDUBNA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliITSsimulation.h"
#include "AliITSresponseSPDdubna.h"
#include "AliITSsegmentationSPD.h"

class AliITSmodule;

//-------------------------------------------------------------------

class AliITSsimulationSPDdubna : public AliITSsimulation {
 public:
    AliITSsimulationSPDdubna();
    AliITSsimulationSPDdubna(AliITSsegmentation *seg, AliITSresponse *res);
    virtual ~AliITSsimulationSPDdubna();
    // copy constructor
    AliITSsimulationSPDdubna(const AliITSsimulationSPDdubna &source); 
     // ass. operator
    AliITSsimulationSPDdubna& operator=(const AliITSsimulationSPDdubna &s);
    // Initilizes the variables
    void Init(AliITSsegmentation *seg, AliITSresponse *resp);

    void InitSimulationModule(Int_t module, Int_t event);
    void SDigitiseModule(AliITSmodule *mod, Int_t mask, Int_t event);
    void WriteSDigits();
    void FinishSDigitiseModule();
    void SDigitsToDigits();
    void DigitiseModule(AliITSmodule *mod,Int_t module,Int_t dummy);
    void UpdateMapSignal(Int_t i,Int_t j,Int_t trk,Int_t ht,Double_t signal);
    void UpdateMapNoise(Int_t ix,Int_t iz,Float_t noise);
    void HitToDigit(AliITSmodule *mod);
    void HitToSDigit(AliITSmodule *mod);
    void HitToSDigitOld(AliITSmodule *mod);
    void ChargeToSignal();
    void CreateHistograms();
    void ResetHistograms();
    TObjArray*  GetHistArray() {return fHis;}// get hist array

 private:
    void SpreadCharge(Double_t x0,Double_t z0,Int_t ix0,Int_t iz0,
		      Double_t el,Double_t sig,Int_t t,Int_t hi);
    AliITSsegmentationSPD* GetSeg(){ // Return pointer to Segmentation class
	return (AliITSsegmentationSPD*)fSegmentation;}
    AliITSresponseSPDdubna* GetResp(){ // Return pointer to Responce class
	return (AliITSresponseSPDdubna*)fResponse;}
    //Double_t * CreateFindCellEdges(Double_t x0,Double_t x1,Double_t z0,
    //                               Double_t z1,Int_t &n);
    void SetCoupling(Int_t row,Int_t col,Int_t ntrack,Int_t idhit);
    void SetCouplingOld(Int_t row, Int_t col,Int_t ntrack,Int_t idhit);
    // Getters for data kept in fSegmentation and fResponse.
    // Returns the Threshold in electrons
    Double_t GetThreshold(){ return GetResp()->GetThreshold();}
    // Returns the couplings Columb and Row.
    void GetCouplings(Float_t &cc,Float_t &cr){
	((AliITSresponseSPDdubna*)GetResp())->GetNoiseParam(cc,cr);}
    // Returns the number of pixels in x
    Int_t GetNPixelsX(){return GetSeg()->Npx();}
    // Returns the number of pixels in z
    Int_t GetNPixelsZ(){return GetSeg()->Npz();}

    ClassDef(AliITSsimulationSPDdubna,2)  // Simulation of SPD clusters

    TObjArray    *fHis;          //! just in case for histogramming

};
#endif 
