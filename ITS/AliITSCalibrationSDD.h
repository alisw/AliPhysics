#ifndef ALIITSCALIBRATIONSDD_H
#define ALIITSCALIBRATIONSDD_H
 
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
/* $Id$ */

#include "AliITSCalibration.h"
#include "AliITSsegmentationSDD.h"
#include "TArrayI.h"

class AliITSCorrMapSDD;
class AliITSDriftSpeedArraySDD;

///////////////////////////////////////////////////////
//  Response for SDD                                 //
///////////////////////////////////////////////////////

class AliITSCalibrationSDD : public AliITSCalibration {
  public:
    //
    // Configuration methods
    //
    AliITSCalibrationSDD();
    AliITSCalibrationSDD(const char *dataType);
    virtual ~AliITSCalibrationSDD();

    virtual Float_t GetBaseline(Int_t anode) const {return fBaseline[anode];}
    virtual void SetBaseline(Int_t anode,Double_t bas) {fBaseline[anode]=bas;}
    virtual Float_t GetNoise(Int_t anode) const {return fNoise[anode];}
    virtual void SetNoise(Int_t anode, Double_t noise) {fNoise[anode]=noise;}

    virtual void  GiveCompressParam(Int_t *x) const;

    void SetZSLowThreshold(Int_t iWing, Int_t thr=25){fZSTL[iWing]=thr;}
    void SetZSHighThreshold(Int_t iWing, Int_t thr=29){fZSTH[iWing]=thr;}
    Int_t GetZSLowThreshold(Int_t iWing) const {return fZSTL[iWing];}
    Int_t GetZSHighThreshold(Int_t iWing) const {return fZSTH[iWing];}

    void  SetNoiseAfterElectronics(Int_t anode,Double_t n=2.38){
	// Noise after electronics (ADC units)
	// 2.36 for ALICE from beam test measurements 2001
	fNoiseAfterEl[anode]=n;}
    Float_t  GetNoiseAfterElectronics(Int_t anode) const {
	// Noise after electronics (ADC units)
	return fNoiseAfterEl[anode];} 
    //void SetDeadChannels(Int_t nchips=0, Int_t nchannels=0);
    void SetDeadChannels(Int_t ndead=0){fDeadChannels=ndead; fBadChannels.Set(ndead);}
    Int_t GetDeadChips() const { return fDeadChips; }
    Int_t GetDeadChannels() const { return fDeadChannels; }
    Float_t GetChannelGain(Int_t anode) const {return fGain[anode];}
    virtual void SetGain(Int_t anode,Double_t g){fGain[anode]=g;}

    
    Int_t GetWing(Int_t anode) const{
      if(anode>=fgkChips*fgkChannels) return 1;
      else return 0;
    }
    Int_t GetChipChannel(Int_t anode) const {return anode%fgkChannels;}
    Int_t GetChip(Int_t anode) const {return anode/fgkChannels;}
    Int_t GetAnodeNumber(Int_t iwing, Int_t ichip03, Int_t ichan) const {
      if(iwing>=2 || ichip03>=4 || ichan>=64) return -1;
      else return iwing*fgkChips*fgkChannels+ichip03*fgkChannels+ichan;
    }
    Int_t GetAnodeNumber(Int_t ichip07, Int_t ichan) const {
      if(ichip07>=8 || ichan>=64) return -1;
      else return ichip07*fgkChannels+ichan;
    }
    
    void    PrintGains() const;
    void    Print();
    virtual void Print(ostream *os) const {AliITSCalibrationSDD::Print(os);}
    virtual void Print(Option_t *option="") const {AliITSCalibrationSDD::Print(option);}
    // not implemented virtual methods (devlared in the mother class
    virtual  void   SetDetParam(Double_t *) 
      {NotImplemented("SetDetParam");}
    virtual void   GetDetParam(Double_t *) const 
      {NotImplemented("GetDetParam");}
    virtual  void   SetNDetParam(Int_t /* n */)
      {NotImplemented("SetNDetParam");}
    virtual Int_t  NDetParam() const
      {NotImplemented("NDetParam"); return 0;}
    virtual void    SetSigmaSpread(Double_t, Double_t) 
      {NotImplemented("SetSigmaSpread");}
    virtual void    SigmaSpread(Double_t & /* p1 */,Double_t & /* p2 */) const 
      {NotImplemented("SigmaSpread");}

    void   SetBad() { 
      fIsBad = kTRUE; 
      for(Int_t i=0;i<fgkChips*fgkWings;i++) fIsChipBad[i]=kTRUE;
    }
    virtual Bool_t IsBad() const { return fIsBad; }
    void   SetChipBad(Int_t nChip) { 
      fIsChipBad[nChip] = kTRUE; 
    }
    virtual Bool_t IsChipBad(Int_t nChip) const { 
      return fIsChipBad[nChip]; 
    }
    Int_t Wings()const{return fgkWings;}//Total number of SDD wings
    Int_t Chips() const{return fgkChips;} // Number of chips/module
    Int_t Channels() const{ return fgkChannels;}//Number of channels/chip
    Int_t NOfAnodes() const {return fgkChannels*fgkChips*fgkWings;}

    virtual void SetBadChannel(Int_t i,Int_t anode);
    Int_t GetBadChannel(Int_t i) const {return fBadChannels[i];}
    Bool_t IsBadChannel(Int_t anode) const{
      if(GetChannelGain(anode)==0) return kTRUE;
      else return kFALSE;
    }
    Float_t GetMapACell(Int_t i,Int_t j) const {
      if(i<256) return fMapAW0->GetCellContent(i,j);
      else return fMapAW1->GetCellContent(i-256,j);
    }
    virtual void SetMapA(Int_t wing,AliITSCorrMapSDD* mapA) {
      if(wing==0) fMapAW0=mapA;
      else fMapAW1=mapA;
    } 
    Float_t GetMapTCell(Int_t i,Int_t j) const {
      if(i<256) return fMapTW0->GetCellContent(i,j);
      else return fMapTW1->GetCellContent(i-256,j);
    }
    virtual void SetMapT(Int_t wing,AliITSCorrMapSDD* mapT) {
      if(wing==0) fMapTW0=mapT;
      else fMapTW1=mapT;
    } 
    
    virtual void SetDriftSpeed(Int_t wing, AliITSDriftSpeedArraySDD* arr){
      if(wing==0) fDrSpeed0=arr;
      else fDrSpeed1=arr;
    }


    virtual Float_t GetDriftSpeedAtAnode(Float_t nAnode) const{
      if(fDrSpeed0==0 || fDrSpeed1==0) AliFatal("Drift speed not set\n");
      if(nAnode<256) return fDrSpeed0->GetDriftSpeed(0,nAnode);
      else return fDrSpeed1->GetDriftSpeed(0,nAnode-256);
    }

    virtual void SetZeroSupp(Bool_t opt=kTRUE) {fZeroSupp=opt;}
    virtual Bool_t GetZeroSupp() const {return fZeroSupp;}

    virtual void SetAMAt40MHz() {fAMAt20MHz=kFALSE;}
    virtual void SetAMAt20MHz() {fAMAt20MHz=kTRUE;}
    virtual Bool_t IsAMAt20MHz() const {return fAMAt20MHz;}

    void GetCorrections(Float_t z, Float_t x, Float_t &devz, Float_t &devx, AliITSsegmentationSDD* seg);
    void GetShiftsForSimulation(Float_t z, Float_t x, Float_t &devz, Float_t &devx, AliITSsegmentationSDD* seg);
    virtual Float_t GetThresholdAnode(Int_t anode, Double_t nsigma=2.2) const {
      return nsigma*fNoiseAfterEl[anode];}


 protected:


    // these statis const should be move to AliITSsegmentationSDD
    static const Int_t fgkWings = 2;     // Number of wings per module
    static const Int_t fgkChips = 4;        // Number of chips/module
    static const Int_t fgkChannels = 64;    // Number of channels/chip
    static const Float_t fgkTemperatureDefault; // default for fT (Kelvin)
    static const Float_t fgkNoiseDefault; // default for fNoise
    static const Float_t fgkBaselineDefault; // default for fBaseline
    static const Float_t fgkGainDefault; //default for gain

    Bool_t fZeroSupp;    // zero suppression
    Bool_t fAMAt20MHz;   // flag for Analog memory of Pascal at 20 MHz
    Int_t fDeadChips;                     // Number of dead chips
    Int_t fDeadChannels;                  // Number of dead channels
    Float_t fGain[fgkWings*fgkChips*fgkChannels];           //Array for channel gains
    Float_t fNoise[fgkWings*fgkChips*fgkChannels];          // Noise array
    Float_t fBaseline[fgkWings*fgkChips*fgkChannels];       // Baseline array
    Float_t fNoiseAfterEl[fgkWings*fgkChips*fgkChannels];   // Noise after electronics

    Int_t   fZSTL[2];     //  Low threshold in 2D zero-suppression (2 hybrids)
    Int_t   fZSTH[2];     // High threshold in 2D zero-suppression (2 hybrids)

    Bool_t   fIsBad;                         // module is dead or alive ?
    Bool_t   fIsChipBad[fgkWings*fgkChips];  // chip is dead or alive ?
    TArrayI  fBadChannels;                   //Array with bad anodes number (0-512) 

    
    AliITSCorrMapSDD* fMapAW0;     //! map of residuals on anode coord. wing 0
    AliITSCorrMapSDD* fMapAW1;     //! map of residuals on anode coord. wing 1
    AliITSCorrMapSDD* fMapTW0;     //! map of residuals on time coord. wing 0
    AliITSCorrMapSDD* fMapTW1;     //! map of residuals on time coord. wing 1
    AliITSDriftSpeedArraySDD* fDrSpeed0; //! drift speed for wing 0
    AliITSDriftSpeedArraySDD* fDrSpeed1; //! drift speed for wing 1

 private:
    AliITSCalibrationSDD(const AliITSCalibrationSDD &ob); // copy constructor
    AliITSCalibrationSDD& operator=(const AliITSCalibrationSDD & /* source */); // ass. op.


    ClassDef(AliITSCalibrationSDD,17) 
    
    };
#endif
