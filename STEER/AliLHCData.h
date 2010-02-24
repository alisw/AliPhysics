#ifndef ALILHCDATA_H
#define ALILHCDATA_H

/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
//  AliLHCData: summary of the lumnosity related information from LHC DIP.     //
//  The time dependent values are averaged over the fPeriod (default: 10 min)  //
//  Created from the TMap provided by the AliLHCReader with optional beginning //
//  and end time stamps to account.                                            //
//                                                                             //
//  The data are (wrapped in the AliLHCDipValT):                               //
//  Total beam intensities   ( eg GetIntensityTotal(beam) for beam=0,1         //
//  Total beam luminosities  ( eg GetLuminTotal(side) for side=0,1 (left/right)//
//  Bunched intensities and luminosities and crossing angle:                   //
//  GetIntensityBunch(beam),GetLuminBunch(side),GetCrossAngle(side)            //
//  Bunches configuration:   GetBunchConfig(beam)                              //
//                                                                             //
//  Collimators information (initial position + changes > 100 microns)         //
//  GetCollimator(collID,jawID) with collID={kTCTVB4L2, kTCTVB4R2, kTCLIA4R2}  //
//  and jawID={kGapDn,kGapUp,kLeftDn,kLeftUp,kRightDn,kRightUp}                //
//                                                                             //
//  Author: ruben.shahoyan@cern.ch                                             //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////
#include "TObject.h"
#include "AliLHCDipValT.h"

class TMap;

class AliLHCData : public TObject
{
 public:
  enum BeamID_t {kBeam1,kBeam2};
  enum Proj_t   {kX,kY};
  enum Side_t   {kLeft,kRight};
  enum Collim_t {kTCTVB4L2, kTCTVB4R2, kTCLIA4R2, kNCollimators};
  enum ColJaw_t {kGapDn,kGapUp,kLeftDn,kLeftUp,kRightDn,kRightUp,kNJaws};
  enum {kRecTotInt,kRecTotIntBunch,kRecBunchInt,kRecTotLum,kRecBunchLum,kRecCrossAngle,kRecBunchConf,kRecFillNum,
	kRecPrfPrID,kRecPrfAmp,kRecPrfPos,kRecPrfSig};
  //
 public:
  //
 AliLHCData() : fPeriod(600),fTMin(0),fTMax(1e20) {}
  AliLHCData(const TMap* dcsMap, double tmin=0, double tmax=1.e20,int avPeriod=600);
  virtual ~AliLHCData() {}
  //
  Bool_t                FillData(const TMap* dcsMap, double tmin=0, double tmax=1.e20);
  //
  const TObjArray&      GetIntensityTotal(BeamID_t b)                const {return fIntTot[b];}
  const TObjArray&      GetIntensityTotalBunch(BeamID_t b)           const {return fIntTotBunch[b];}
  const TObjArray&      GetIntensityBunch(BeamID_t b)                const {return fIntBunch[b];}
  const TObjArray&      GetBeamPos(BeamID_t b,Proj_t p)              const {return fBeamPos[b][p];}
  const TObjArray&      GetLuminTotal(Side_t s)                      const {return fLuminTot[s];}
  const TObjArray&      GetLuminBunch(Side_t s)                      const {return fLuminBunch[s];}
  const TObjArray&      GetCrossAngle(Side_t s)                      const {return fCrossAngle[s];}
  const AliLHCDipValI&  GetBunchConfig(BeamID_t b)                   const {return fBunchConfig[b];}
  const TObjArray&      GetCollimator(Collim_t c,ColJaw_t jaw)       const {return fCollimators[c][jaw];}
  //
  Double_t              GetTMin()                                    const {return fTMin;}
  Double_t              GetTMax()                                    const {return fTMax;}
  Int_t                 GetPeriod()                                  const {return fPeriod;}
  Int_t                 GetNBunches(BeamID_t b)                      const {return GetBunchConfig(b).GetSize();}
  Int_t                 GetBunchRFBeam1(BeamID_t b,Int_t bunch)      const {return GetBunchConfig(b).GetValue(bunch);}
  Int_t                 GetFillNumber()                              const {return GetUniqueID();}
  void                  SetFillNumber(Int_t fill)                          {SetUniqueID(fill);}
  void                  SetPeriod(Int_t nsec=600)                          {fPeriod = nsec;}
  void                  SetTMin(Double_t t)                                {fTMin = t;}
  void                  SetTMax(Double_t t)                                {fTMax = t;}
  virtual void          Print(const Option_t *opt="")                const;
  //
 protected:
  TObjArray*            GetDCSEntry(const TMap* dcsMap,const char* key,int &entry,double tmin,double tmax) const;

 AliLHCData(const AliLHCData& src) : TObject(src),fPeriod(src.fPeriod),fTMin(src.fTMin),fTMax(src.fTMax) { /*dummy*/ }
  AliLHCData& operator=(const AliLHCData& ) { /*dummy*/ return *this;}
  
 protected:
  //
  Int_t           fPeriod;                            // averaging period in seconds
  Double_t        fTMin;                              // selection timeMin
  Double_t        fTMax;                              // selection timeMax
  AliLHCDipValI   fBunchConfig[2];                    // bunches configuration for beam1 and beam2 
  //
  TObjArray       fIntTot[2];                         // total intensity (DC BCT) for beam1 and beam2 (AliLHCDipValD)
  TObjArray       fIntTotBunch[2];                    // total intensity (BCTFR) for beam1 and beam2 (AliLHCDipValD)
  TObjArray       fIntBunch[2];                       // bunch by bunch intensity for beam1 and beam2 (AliLHCDipValD)
  //
  TObjArray       fBeamPos[2][2];                     // horizontal and vertical projection gaussian fit params for beam1 and beam2 (AliLHCDipValD: amp,pos,sigma)
  //
  TObjArray       fLuminTot[2];                       // total luminosity from BRANB_4L2 and BRANB_4R2 (AliLHCDipValD)
  TObjArray       fLuminBunch[2];                     // bunch by bunch luminosity from BRANB_4L2 and BRANB_4R2 (AliLHCDipValD)
  TObjArray       fCrossAngle[2];                     // beams crossing angle from BRANB_4L2 and BRANB_4R2 (AliLHCDipValD)
  // 
  TObjArray       fCollimators[kNCollimators][kNJaws];// collimators data (AliLHCDipValD: kGapDn,kGapUp,kLeftDn,kLeftUp,kRightDn,kRightUp)
  //
  static const Char_t *fgkDCSNames[];                 // beam related DCS names to extract
  static const Char_t *fgkDCSColNames[];              // collimators to extract
  static const Char_t *fgkDCSColJaws[];               // names of collimator pieces

  ClassDef(AliLHCData,2)
};

#endif
