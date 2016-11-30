#ifndef ALILHCDATA_H
#define ALILHCDATA_H

/********************************************************************************
*                                                                               *
*   AliLHCData: summary of the LHC related information from LHC DIP.            *
*   Created from the TMap provided by the AliLHCReader with optional beginning  *
*                                                                               *
*   The data are (wrapped in the AliLHCDipValT):                                *
*   made of TimeStamp (double) and array of values                              *
*                                                                               *
*   Multiple entries for each type of data are possible. To obtaine number of   *
*   records (with distinct timestamp) for give type od records use:             *
*   int GetNBunchConfigMeasured(int beam) (with beam=0,1) etc.                  *
*                                                                               *
*   To get i-th entry, use brec= AliLHCDipValI* GetBunchConfigMeasured(bm,i);   *
*   Note: exact type of templated AliLHCDipValT pointer depends on the record   *
*   type, concult getters to know it.                                           *
*                                                                               *
*   Then, once the pointer is obtained, details can be accessed:                *
*   int nBunches = brec->GetSize();                                             *
*   for (int i=0;i<nBunches;i++) printf("Bunch#%d: %d\n",i,(*brec)[i]);         *
*                                                                               *
*   ATTENTION: Bunch RFBucked is NEGATIVE for bunches interacting at IR2        *
*                                                                               *
*                                                                               *
*                                                                               *
*   Author: ruben.shahoyan@cern.ch                                              *
*                                                                               *
********************************************************************************/

#include "AliLHCDipValT.h"
#include "TObject.h"
class TObjArray;
//class AliLHCDipValT;

class AliDCSArray;
class TString;
class TMap;
class AliLHCReader;
class TGraph;

class AliLHCData : public TObject
{
 public:
  enum          {kStart,kNStor};
  enum BeamID_t {kBeam1,kBeam2};
  enum BgID_t   {kBg1,kBg2,kBg3,kNBGs};
  enum Proj_t   {kX,kY};
  enum Side_t   {kLeft,kRight};
  enum Collim_t {kTCTVB4L2, kTCTVB4R2, kTCLIA4R2, kNCollimators};
  enum ColJaw_t {kGapDn,kGapUp,kLeftDn,kLeftUp,kRightDn,kRightUp,kNJaws};
  enum          {kMaxBSlots = 3564, kOffsBeam1=346, kOffsBeam2 = 3019, kY2015 = 1420070400};
  enum          {kMarginSOR = 60*60*24*30, // use margin of 30 days for SOR, when looking for the 1st record
		 kMarginEOR = 60*15};      // use margin of 15 min for EOR, when looking for the last record
  //
  enum {kIntTot,kIntTotAv,kIntBunchAv,
	kLumAcqMode,kLumTot,kLumTotErr,kLumBunch,kLumBunchErr,kLumCrossAng,kLumCrossAngErr,
	kBunchConf,kFillNum,kBunchLgtNB,kBunchLgt,kBunchLgtFillB,
	kRCInjSch,kRCBeta,kRCCrossAng,kRCVang,
	kBeamSzAcqMode,kBeamSzSigH,kBeamSzSigV,kBeamSzEmittH,kBeamSzEmittV,kBeamSzSigHErr,kBeamSzSigVErr,
	kCollPos,
	kBPTXdeltaTB1B2,
	kBPTXdeltaTRMSB1B2,
	kBPTXPhase,
	kBPTXPhaseRMS,
	kBPTXPhaseShift,
	//
	kALILumiTotalInst,
	kALILumiTotalDeliveredStabBeam,
	kALILumiBunchInst,
	kALIBackground,
	kNRecordTypes};
  //
  //le
 public:
  //
 AliLHCData() : fTMin(0),fTMax(1e10),fFillNumber(0),fData(0),fkFile2Process(0),fkMap2Process(0) {Clear();}
  AliLHCData(const TMap*   dcsMap,  double tmin=0, double tmax=1.e10);
  AliLHCData(const Char_t* dcsFile, double tmin=0, double tmax=1.e10);
  virtual ~AliLHCData() {}
  //
  Bool_t                FillData(const TMap*   dcsMap,  double tmin=0, double tmax=1.e20);
  Bool_t                FillData(const Char_t* dcsFile, double tmin=0, double tmax=1.e20);
  Double_t              GetTMin()                                    const {return fTMin;}
  Double_t              GetTMax()                                    const {return fTMax;}
  Int_t                 GetFillNumber()                              const {return fFillNumber;}
  void                  SetFillNumber(Int_t fill)                          {fFillNumber = fill;}
  void                  SetTMin(Double_t t)                                {fTMin = t<0?0:(t>1e10?1e10:t);}
  void                  SetTMax(Double_t t)                                {fTMax = t<0?0:(t>1e10?1e10:t);}
  Bool_t                IsRun2()                                     const {return fTMin>kY2015;}
  //
  virtual void          Print(const Option_t *opt="")                const;
  TGraph*               ExportGraph(Int_t *coord, Int_t elID=0)      const;
  //
  Int_t GetNBunchConfigMeasured(int bm)           const {return GoodPairID(bm)?fBunchConfMeas[bm][kNStor]:-1;}
  Int_t GetNBunchConfigDeclared(int bm)           const {return GoodPairID(bm)?fBunchConfDecl[bm][kNStor]:-1;}
  Int_t GetNBunchLengths(int bm)                  const {return GoodPairID(bm)?fBunchLengths[bm][kNStor]:-1;}
  Int_t GetNTotalIntensity(int bm)                const {return GoodPairID(bm)?fIntensTotal[bm][kNStor]:-1;}
  Int_t GetNTotalIntensityAv(int bm)              const {return GoodPairID(bm)?fIntensTotalAv[bm][kNStor]:-1;}
  Int_t GetNIntensityPerBunch(int bm)             const {return GoodPairID(bm)?fIntensPerBunch[bm][kNStor]:-1;}
  Int_t GetNEmittanceH(int bm)                    const {return GoodPairID(bm)?fEmittanceH[bm][kNStor]:-1;}
  Int_t GetNEmittanceV(int bm)                    const {return GoodPairID(bm)?fEmittanceV[bm][kNStor]:-1;}
  Int_t GetNBeamSigmaH(int bm)                    const {return GoodPairID(bm)?fBeamSigmaH[bm][kNStor]:-1;}
  Int_t GetNBeamSigmaV(int bm)                    const {return GoodPairID(bm)?fBeamSigmaV[bm][kNStor]:-1;}
  //
  Int_t GetNLuminosityTotal(int lr)               const {return GoodPairID(lr)?fLuminTotal[lr][kNStor]:-1;}
  Int_t GetNLuminosityPerBunch(int lr)            const {return GoodPairID(lr)?fLuminPerBC[lr][kNStor]:-1;}
  Int_t GetNLuminosityAcqMode(int lr)             const {return GoodPairID(lr)?fLuminAcqMode[lr][kNStor]:-1;}
  Int_t GetNCrossingAngle(int lr)                 const {return GoodPairID(lr)?fCrossAngle[lr][kNStor]:-1;}
  //
  Int_t GetNInjectionScheme()                     const {return fRCInjScheme[kNStor];}
  Int_t GetNRCBetaStar()                          const {return fRCBeta[kNStor];}
  Int_t GetNRCAngleH()                            const {return fRCAngH[kNStor];}
  Int_t GetNRCAngleV()                            const {return fRCAngV[kNStor];}
  //
  Int_t GetNLumiAlice()                           const {return fLumiAlice[kNStor];}
  Int_t GetNLumiAliceSBDelivered()                const {return fLumiAliceStB[kNStor];}
  Int_t GetNLumiAliceBunch()                      const {return fLumiAliceBbB[kNStor];}  
  //
  Int_t GetNBckgAlice(int bg)                     const {return (bg>-1&&bg<kNBGs) ? fBckgAlice[bg][kNStor] : -1;}
  //
  Int_t GetNBPTXdeltaTB1B2()                      const {return fBPTXdTB1B2[kNStor];}
  Int_t GetNBPTXdeltaTRMSB1B2()                   const {return fBPTXdTRMSB1B2[kNStor];}
  Int_t GetNBPTXPhase(int bm)                     const {return GoodPairID(bm)?fBPTXPhase[bm][kNStor]:-1;}
  Int_t GetNBPTXPhaseRMS(int bm)                  const {return GoodPairID(bm)?fBPTXPhaseRMS[bm][kNStor]:-1;}
  Int_t GetNBPTXPhaseShift(int bm)                const {return GoodPairID(bm)?fBPTXPhaseShift[bm][kNStor]:-1;}
  //
  //
  Int_t GetNCollimatorJawPos(int coll,int jaw)    const;
  //
  AliLHCDipValI* GetBunchConfigMeasured(int bm, int i=0)  const;
  AliLHCDipValF* GetBunchLengths(int bm, int i=0)         const;
  AliLHCDipValI* GetBunchConfigDeclared(int bm, int i=0)  const;
  AliLHCDipValF* GetTotalIntensity(int bm, int i=0)       const;
  AliLHCDipValF* GetTotalIntensityAv(int bm, int i=0)     const;
  AliLHCDipValF* GetIntensityPerBunch(int bm, int i=0)    const;
  AliLHCDipValF* GetEmittanceH(int bm, int i=0)           const;
  AliLHCDipValF* GetEmittanceV(int bm, int i=0)           const;
  AliLHCDipValF* GetBeamSigmaH(int bm, int i=0)           const;
  AliLHCDipValF* GetBeamSigmaV(int bm, int i=0)           const;
  AliLHCDipValF* GetLuminosityTotal(int lr, int i=0)      const;
  AliLHCDipValF* GetLuminosityPerBunch(int lr, int i=0)   const;
  AliLHCDipValI* GetLuminosityAcqMode(int lr, int i=0)    const;
  AliLHCDipValF* GetCrossAngle(int lr, int i=0)           const;
  AliLHCDipValC* GetInjectionScheme(int i=0)              const;
  AliLHCDipValF* GetRCBetaStar(int i=0)                   const;
  AliLHCDipValF* GetRCAngleH(int i=0)                     const; 
  AliLHCDipValF* GetRCAngleV(int i=0)                     const; 
  AliLHCDipValF* GetCollimJawPos(int coll, int jaw, int i=0) const;
  //
  AliLHCDipValF* GetLumiAlice(int i=0)                    const; 
  AliLHCDipValF* GetLumiAliceSBDelivered(int i=0)         const; 
  AliLHCDipValF* GetLumiAliceBunch(int i=0)               const; 
  AliLHCDipValF* GetBckgAlice(int bg, int i=0)            const; 
  //
  AliLHCDipValF* GetBPTXdeltaTB1B2(int i=0)               const;
  AliLHCDipValF* GetBPTXdeltaTRMSB1B2(int i=0)            const;
  AliLHCDipValF* GetBPTXPhase(int bm, int i=0)            const;
  AliLHCDipValF* GetBPTXPhaseRMS(int bm, int i=0)         const;
  AliLHCDipValF* GetBPTXPhaseShift(int bm, int i=0)       const;
  //
  Float_t        GetLumiAlice(Double_t tstamp)            const;
  Float_t        GetLumiAliceSBDelivered(Double_t tstamp) const;
  Float_t        GetBckgAlice(Int_t bg,Double_t tstamp)   const;
  //
  //  Float_t        GetLumiInstAlice(Double_t tstamp)        const;
  //  Float_t        GetBckgInstAlice(Double_t tstamp)        const;
  //
  void           FlagInteractingBunches(const Int_t beam1[2],const Int_t beam2[2]);
  TObject*       FindRecValidFor(int start,int nrec, double tstamp) const;
  Int_t          FindEntryValidFor(int start,int nrec, double tstamp) const;
  AliLHCDipValI* GetBunchConfigMeasured(int beam,double tstamp)  const;
  AliLHCDipValI* GetBunchConfigDeclared(int beam,double tstamp)  const;
  Int_t          GetNInteractingBunchesMeasured(int i=0)         const;
  Int_t          GetNInteractingBunchesDeclared(int i=0)         const;
  Int_t          IsPilotPresent(int i=0)                         const;
  //
  // return array with beginning [0] and number of records for corresponding info (in the fData)
  const Int_t* GetOffsBunchConfigMeasured(int bm)         const {return GoodPairID(bm)?fBunchConfMeas[bm]:0;}
  const Int_t* GetOffsBunchConfigDeclared(int bm)         const {return GoodPairID(bm)?fBunchConfDecl[bm]:0;}
  const Int_t* GetOffsBunchLengths(int bm)                const {return GoodPairID(bm)?fBunchLengths[bm]:0;}
  const Int_t* GetOffsTotalIntensity(int bm)              const {return GoodPairID(bm)?fIntensTotal[bm]:0;}
  const Int_t* GetOffsTotalIntensityAv(int bm)            const {return GoodPairID(bm)?fIntensTotalAv[bm]:0;}
  const Int_t* GetOffsIntensityPerBunch(int bm)           const {return GoodPairID(bm)?fIntensPerBunch[bm]:0;}
  const Int_t* GetOffsEmittanceH(int bm)                  const {return GoodPairID(bm)?fEmittanceH[bm]:0;}
  const Int_t* GetOffsEmittanceV(int bm)                  const {return GoodPairID(bm)?fEmittanceV[bm]:0;}
  const Int_t* GetOffsBeamSigmaH(int bm)                  const {return GoodPairID(bm)?fBeamSigmaH[bm]:0;}
  const Int_t* GetOffsBeamSigmaV(int bm)                  const {return GoodPairID(bm)?fBeamSigmaV[bm]:0;}
  //
  const Int_t* GetOffsLuminosityTotal(int lr)             const {return GoodPairID(lr)?fLuminTotal[lr]:0;}
  const Int_t* GetOffsLuminosityPerBunch(int lr)          const {return GoodPairID(lr)?fLuminPerBC[lr]:0;}
  const Int_t* GetOffsLuminosityAcqMode(int lr)           const {return GoodPairID(lr)?fLuminAcqMode[lr]:0;}
  const Int_t* GetOffsCrossingAngle(int lr)               const {return GoodPairID(lr)?fCrossAngle[lr]:0;}
  //
  const Int_t* GetOffsInjectionScheme()                   const {return fRCInjScheme;}
  const Int_t* GetOffsRCBetaStar()                        const {return fRCBeta;}
  const Int_t* GetOffsRCAngleH()                          const {return fRCAngH;}
  const Int_t* GetOffsRCAngleV()                          const {return fRCAngV;}
  //
  const Int_t* GetOffsLumiAlice()                         const {return fLumiAlice;}
  const Int_t* GetOffsLumiAliceSBDelivered()              const {return fLumiAliceStB;}
  const Int_t* GetOffsLumiAliceBunch()                    const {return fLumiAliceBbB;}
  const Int_t* GetOffsBckgAlice(int bg)                   const {return (bg>-1&&bg<kNBGs) ? fBckgAlice[bg] : 0;}
  //
  const Int_t* GetOffsBPTXdeltaTB1B2()                       const {return fBPTXdTB1B2;}
  const Int_t* GetOffsBPTXdeltaTRMSB1B2()                    const {return fBPTXdTRMSB1B2;}
  const Int_t* GetOffsBPTXPhase(int bm)                      const {return GoodPairID(bm)?fBPTXPhase[bm]:0;}
  const Int_t* GetOffsBPTXPhaseRMS(int bm)                   const {return GoodPairID(bm)?fBPTXPhaseRMS[bm]:0;}
  const Int_t* GetOffsBPTXPhaseShift(int bm)                 const {return GoodPairID(bm)?fBPTXPhaseShift[bm]:0;}
  //
  const Int_t* GetOffsCollimatorJawPos(int coll,int jaw)  const;
  //
  const TObjArray&  GetData()                             const {return fData;}
  //
  // analysis methods
  Int_t GetMeanIntensity(int beamID, Double_t &colliding, Double_t &noncolliding, const TObjArray* bcmasks=0) const;
  Int_t GetBCId(int bucket, int beamID) const;
  //
  // for retrofitting, these methods has to be public
  //  void                  FillLumiAliceOFL(Int_t nrec, Int_t* time, Double_t* val);
  //  void                  FillBckgAliceOFL(Int_t nrec, Int_t* time, Double_t* val);

 protected:
  //
  Bool_t                FillData(double tmin=0, double tmax=1.e20);
  virtual void          Clear(const Option_t *opt="");
  void                  PrintAux(Bool_t full,const Int_t refs[2],const Option_t *opt="") const;
  TObjArray*            GetDCSEntry(const char* key,int &entry,int &last,double tmin,double tmax) const;
  Int_t                 FillScalarRecord(  int refs[2], const char* rec, const char* recErr=0, Double_t maxAbsVal=1.e30);
  Int_t                 FillBunchConfig(   int refs[2], const char* rec);
  Int_t                 FillStringRecord(  int refs[2], const char* rec);
  Int_t                 FillAcqMode(       int refs[2], const char* rec);
  Int_t                 FillBunchInfo(     int refs[2], const char* rec,int ibm, Bool_t inRealSlots, Double_t maxAbsVal=1.e30);
  Int_t                 FillBCLuminosities(int refs[2], const char* rec, const char* recErr, Int_t useBeam, Double_t maxAbsVal=1.e30);
  //
  Int_t                 ExtractInt(AliDCSArray* dcsArray,Int_t el)    const;
  Double_t              ExtractDouble(AliDCSArray* dcsArray,Int_t el) const;
  TString&              ExtractString(AliDCSArray* dcsArray)          const;
  AliLHCData(const AliLHCData& src) : TObject(src),fTMin(0),fTMax(0),fFillNumber(0),fData(0),fkFile2Process(0),fkMap2Process(0) { /*dummy*/ }
  AliLHCData& operator=(const AliLHCData& ) { /*dummy*/ return *this;}
  Int_t                 TimeDifference(double v1,double v2,double tol=0.9) const;
  Bool_t                IzZero(double val, double tol=1e-16)         const {return TMath::Abs(val)<tol;}
  Bool_t                GoodPairID(int beam)                         const;
  //
 protected:
  //
  Double_t        fTMin;                              // selection timeMin
  Double_t        fTMax;                              // selection timeMax
  Int_t           fFillNumber;                        // fill number           : kFillNum
  //
  //---------------- Last index gives: 0 - beginning of the records in fData, 1 - number of records
  //
  //              infrormation from RunControl
  Int_t           fRCInjScheme[2];                    // active injection scheme                       : String |kRCInjScheme
  Int_t           fRCBeta[2];                         // target beta                                   : Float  |kRCBeta
  Int_t           fRCAngH[2];                         // horisontal angle                              : Float  |kRCCrossAng
  Int_t           fRCAngV[2];                         // vertical angle                                : Float  |kRCVang
  Int_t           fBunchConfDecl[2][2];               // declared beam configuration                   : Float  |kBunchConf                
  //
  //              measured information
  Int_t           fBunchConfMeas[2][2];               // measured beam configuration                   : Int    |kBunchLgtFillB
  Int_t           fBunchLengths[2][2];                // measured beam lenghts                         : Float  |kBunchLgt
  Int_t           fIntensTotal[2][2];                 // total beam intensities                        : Float  |kIntTot
  Int_t           fIntensTotalAv[2][2];               // total beam intensities from bunch averages    : Float  |kIntTotAv
  Int_t           fIntensPerBunch[2][2];              // bunch-by-bunch intensities                    : Float  |kIntBunchAv
  //
  Int_t           fCrossAngle[2][2];                  // crossing angle   at IP2 and its error         : Float  |kLimCrossAng, kLumCrossAngErr
  Int_t           fEmittanceH[2][2];                  // beam H emittances                             : Float  |kBeamSzEmittH
  Int_t           fEmittanceV[2][2];                  // beam V emittances                             : Float  |kBeamSzEmittV
  Int_t           fBeamSigmaH[2][2];                  // beam H sigma and error                        : Float  |kBeamSzSigH,kBeamSzSigHErr
  Int_t           fBeamSigmaV[2][2];                  // beam V sigma and error                        : Float  |kBeamSzSigV,kBeamSzSigVErr
  //
  Int_t           fLuminTotal[2][2];                  // total luminosity at IP2 and its error         : Float  |kLumTot, kLumTotErr
  Int_t           fLuminPerBC[2][2];                  // luminosity at IP2 for each BC and its error   : Float  |kLumBunch,kLumBunchErr
  Int_t           fLuminAcqMode[2][2];                // luminosity acquisition mode                   : Int    |kLumAcqMode
  //
  // here we will store the luminosity and the background measured by Alice. We store the value integrated from the start of fill.
  // the inst. value can be obtained as its derivative
  Int_t           fLumiAlice[2];                      // luminosity measured by Alice, Total_Inst      : Float  |kALILumiTotalInst
  Int_t           fLumiAliceStB[2];                   // luminosity measured by Alice, Deli.StabBeam   : Float  |kALILumiTotalDeliveredStabBeam
  Int_t           fLumiAliceBbB[2];                   // luminosity measured by Alice, B-by-B          : Float  |kALILumiBunchInst
  //
  Int_t           fBckgAlice[kNBGs][2];               // backgrounds measured by Alice                 : Float  |kALIBackground1,2,3
  //
  Int_t           fBPTXdTB1B2[2];                     // BPTX_deltaT_B1_B2                             : Float  |kBPTXdeltaTB1B2,
  Int_t           fBPTXdTRMSB1B2[2];                  // BPTX_deltaTRMS_B1_B2                          : Float  |kBPTXdeltaTRMSB1B2,
  Int_t           fBPTXPhase[2][2];                   // BPTX_Phase_B1, and BPTX_Phase_B2              : Float  |kBPTXPhase
  Int_t           fBPTXPhaseRMS[2][2];                // BPTX_PhaseRMS_B1, and BPTX_PhaseRMS_B2        : Float  |kBPTXPhaseRMS
  Int_t           fBPTXPhaseShift[2][2];              // BPTX_Phase_Shift_B1 and BPTX_Phase_Shift_B2   : Float  |kBPTXPhaseShift 
  //
  Int_t           fCollimators[kNCollimators][kNJaws][2];// collimator jaws positions                  : Float  |kCollPos
  //
  TObjArray       fData;                              // single storage for various records
  //
  static const Char_t *fgkDCSNames[];                 // beam related DCS names to extract
  static const Char_t *fgkDCSColNames[];              // collimators to extract
  static const Char_t *fgkDCSColJaws[];               // names of collimator pieces
  //
 private:
  // non-persistent objects used at the filling time
  const Char_t*   fkFile2Process;                      //! name of DCS file
  const TMap*     fkMap2Process;                       //! DCS map to process 

  ClassDef(AliLHCData,4)
};


//_____________________________________________________________________________
inline Int_t AliLHCData::GetNCollimatorJawPos(int coll,int jaw) const {// get n records
  return (coll>=0&&coll<kNCollimators&&jaw>=0&&jaw<kNJaws)? fCollimators[coll][jaw][kNStor]:0;
}

inline const Int_t* AliLHCData::GetOffsCollimatorJawPos(int coll,int jaw)  const { // offset array
  return (coll>=0&&coll<kNCollimators&&jaw>=0&&jaw<kNJaws)? fCollimators[coll][jaw]:0;
}

inline AliLHCDipValI* AliLHCData::GetBunchConfigMeasured(int bm, int i) const { // get record
  return (GoodPairID(bm) && i>=0 && i<fBunchConfMeas[bm][kNStor]) ? (AliLHCDipValI*)fData[fBunchConfMeas[bm][kStart]+i]:0;
}

inline AliLHCDipValF* AliLHCData::GetBunchLengths(int bm, int i) const { // get record
  return (GoodPairID(bm) && i>=0 && i<fBunchLengths[bm][kNStor]) ? (AliLHCDipValF*)fData[fBunchLengths[bm][kStart]+i]:0;
}

inline AliLHCDipValI* AliLHCData::GetBunchConfigDeclared(int bm, int i) const { // get record
  return (GoodPairID(bm) && i>=0 && i<fBunchConfDecl[bm][kNStor]) ? (AliLHCDipValI*)fData[fBunchConfDecl[bm][kStart]+i]:0;
}

inline AliLHCDipValF* AliLHCData::GetTotalIntensity(int bm, int i) const { // get record
  return (GoodPairID(bm) && i>=0 && i<fIntensTotal[bm][kNStor]) ? (AliLHCDipValF*)fData[fIntensTotal[bm][kStart]+i]:0;
}

inline AliLHCDipValF* AliLHCData::GetTotalIntensityAv(int bm, int i) const { // get record
  return (GoodPairID(bm) && i>=0 && i<fIntensTotalAv[bm][kNStor]) ? (AliLHCDipValF*)fData[fIntensTotalAv[bm][kStart]+i]:0;
}

inline AliLHCDipValF* AliLHCData::GetIntensityPerBunch(int bm, int i) const { // get record
  return (GoodPairID(bm) && i>=0 && i<fIntensPerBunch[bm][kNStor]) ? (AliLHCDipValF*)fData[fIntensPerBunch[bm][kStart]+i]:0;
}

inline AliLHCDipValF* AliLHCData::GetEmittanceH(int bm, int i) const { // get record
  return (GoodPairID(bm) && i>=0 && i<fEmittanceH[bm][kNStor]) ? (AliLHCDipValF*)fData[fEmittanceH[bm][kStart]+i]:0;
}

inline AliLHCDipValF* AliLHCData::GetEmittanceV(int bm, int i) const { // get record
  return (GoodPairID(bm) && i>=0 && i<fEmittanceV[bm][kNStor]) ? (AliLHCDipValF*)fData[fEmittanceV[bm][kStart]+i]:0;
}

inline AliLHCDipValF* AliLHCData::GetBeamSigmaH(int bm, int i) const { // get record
  return (GoodPairID(bm) && i>=0 && i<fBeamSigmaH[bm][kNStor]) ? (AliLHCDipValF*)fData[fBeamSigmaH[bm][kStart]+i]:0;
}

inline AliLHCDipValF* AliLHCData::GetBeamSigmaV(int bm, int i) const { // get record
  return (GoodPairID(bm) && i>=0 && i<fBeamSigmaV[bm][kNStor]) ? (AliLHCDipValF*)fData[fBeamSigmaV[bm][kStart]+i]:0;
}

inline AliLHCDipValF* AliLHCData::GetLuminosityTotal(int lr, int i) const { // get record
  return (GoodPairID(lr) && i>=0 && i<fLuminTotal[lr][kNStor]) ? (AliLHCDipValF*)fData[fLuminTotal[lr][kStart]+i]:0;
}

inline AliLHCDipValF* AliLHCData::GetLuminosityPerBunch(int lr, int i) const { // get record
  return (GoodPairID(lr) && i>=0 && i<fLuminPerBC[lr][kNStor]) ? (AliLHCDipValF*)fData[fLuminPerBC[lr][kStart]+i]:0;
}

inline AliLHCDipValI* AliLHCData::GetLuminosityAcqMode(int lr, int i) const { // get record
  return (GoodPairID(lr) && i>=0 && i<fLuminAcqMode[lr][kNStor]) ? (AliLHCDipValI*)fData[fLuminAcqMode[lr][kStart]+i]:0;
}

inline AliLHCDipValF* AliLHCData::GetCrossAngle(int lr, int i) const { // get record
  return (GoodPairID(lr) && i>=0 && i<fCrossAngle[lr][kNStor]) ? (AliLHCDipValF*)fData[fCrossAngle[lr][kStart]+i]:0;
}

inline AliLHCDipValC* AliLHCData::GetInjectionScheme(int i) const { // get record
  return (i>=0 && i<fRCInjScheme[kNStor]) ? (AliLHCDipValC*)fData[fRCInjScheme[kStart]+i]:0;
}

inline AliLHCDipValF* AliLHCData::GetRCBetaStar(int i) const { // get record
  return (i>=0 && i<fRCBeta[kNStor]) ? (AliLHCDipValF*)fData[fRCBeta[kStart]+i]:0;
}

inline AliLHCDipValF* AliLHCData::GetRCAngleH(int i) const { // get record
  return (i>=0 && i<fRCAngH[kNStor]) ? (AliLHCDipValF*)fData[fRCAngH[kStart]+i]:0;
}

inline AliLHCDipValF* AliLHCData::GetRCAngleV(int i) const { // get record
  return (i>=0 && i<fRCAngV[kNStor]) ? (AliLHCDipValF*)fData[fRCAngV[kStart]+i]:0;
}

inline AliLHCDipValF* AliLHCData::GetCollimJawPos(int coll, int jaw, int i) const { // get record
  return (coll>=0 && coll<kNCollimators && jaw>=0 && jaw<kNJaws && 
	  i>=0 && i<fCollimators[coll][jaw][kNStor]) ? (AliLHCDipValF*)fData[fCollimators[coll][jaw][kStart]+i]:0;
}

inline AliLHCDipValF* AliLHCData::GetLumiAlice(int i) const { // get record on integrated luminosity
  return (i>=0 && i<fLumiAlice[kNStor]) ? (AliLHCDipValF*)fData[fLumiAlice[kStart]+i]:0;
}

inline AliLHCDipValF* AliLHCData::GetLumiAliceSBDelivered(int i) const { // get record on st.beam delivered luminosity
  return (i>=0 && i<fLumiAliceStB[kNStor]) ? (AliLHCDipValF*)fData[fLumiAliceStB[kStart]+i]:0;
}

inline AliLHCDipValF* AliLHCData::GetLumiAliceBunch(int i) const { // get record on b-by-b luminosity
  return (i>=0 && i<fLumiAliceBbB[kNStor]) ? (AliLHCDipValF*)fData[fLumiAliceBbB[kStart]+i]:0;
}

inline AliLHCDipValF* AliLHCData::GetBckgAlice(int bg, int i) const { // get record on integrated background
  return (bg>=0&&bg<kNBGs&&i>=0&&i<fBckgAlice[bg][kNStor]) ? (AliLHCDipValF*)fData[fBckgAlice[bg][kStart]+i]:0;
}

inline AliLHCDipValF* AliLHCData::GetBPTXdeltaTB1B2(int i) const { // get record 
  return (i>=0 && i<fBPTXdTB1B2[kNStor]) ? (AliLHCDipValF*)fData[fBPTXdTB1B2[kStart]+i]:0;
}

inline AliLHCDipValF* AliLHCData::GetBPTXdeltaTRMSB1B2(int i) const { // get record 
  return (i>=0 && i<fBPTXdTRMSB1B2[kNStor]) ? (AliLHCDipValF*)fData[fBPTXdTRMSB1B2[kStart]+i]:0;
}

inline AliLHCDipValF* AliLHCData::GetBPTXPhase(int bm, int i) const { // get record 
  return (GoodPairID(bm) && i>=0 && i<fBPTXPhase[bm][kNStor]) ? (AliLHCDipValF*)fData[fBPTXPhase[bm][kStart]+i]:0;
}

inline AliLHCDipValF* AliLHCData::GetBPTXPhaseRMS(int bm, int i) const { // get record 
  return (GoodPairID(bm) && i>=0 && i<fBPTXPhaseRMS[bm][kNStor]) ? (AliLHCDipValF*)fData[fBPTXPhaseRMS[bm][kStart]+i]:0;
}

inline AliLHCDipValF* AliLHCData::GetBPTXPhaseShift(int bm, int i) const { // get record 
  return (GoodPairID(bm) && i>=0 && i<fBPTXPhaseShift[bm][kNStor]) ? (AliLHCDipValF*)fData[fBPTXPhaseShift[bm][kStart]+i]:0;
}

inline Float_t AliLHCData::GetLumiAlice(Double_t tStamp) const { // get closest in time value on integrated luminosity
  int idx = FindEntryValidFor(fLumiAlice[kStart],fLumiAlice[kNStor],tStamp);
  return idx<0 ? -1 : ((AliLHCDipValF*)fData[fLumiAlice[kStart]+idx])->GetValue();
}

inline Float_t AliLHCData::GetLumiAliceSBDelivered(Double_t tStamp) const { // get closest in time value on delivered luminosity
  int idx = FindEntryValidFor(fLumiAliceStB[kStart],fLumiAliceStB[kNStor],tStamp);
  return idx<0 ? -1 : ((AliLHCDipValF*)fData[fLumiAliceStB[kStart]+idx])->GetValue();
}

inline Float_t AliLHCData::GetBckgAlice(int bg,Double_t tStamp) const { // get closest in time value on integrated bckg
  if (bg<0||bg>=kNBGs) return 0;
  int idx = FindEntryValidFor(fBckgAlice[bg][kStart],fBckgAlice[bg][kNStor],tStamp);
  return idx<0 ? -1 : ((AliLHCDipValF*)fData[fBckgAlice[bg][kStart]+idx])->GetValue();
}



inline Int_t AliLHCData::FindEntryValidFor(int start,int nrec, double tstamp) const
{
  // find index of record within this limits valid for given tstamp (i.e. the last one before or equal to tstamp)
  int idx;
  for (idx=0;idx<nrec;idx++) {
    if (TimeDifference(tstamp,((AliLHCDipValI*)fData[start+idx])->GetTimeStamp())<=0) break;
  }
  return (idx<nrec) ? idx : nrec-1;
}


inline Int_t AliLHCData::GetBCId(int bucket, int beamID) const
{
  Int_t offs = beamID==0 ? kOffsBeam1:kOffsBeam2;
  if (IsRun2()) offs -= 2; // constants were changed (Martino)
  return (TMath::Abs(bucket)/10 + offs)%kMaxBSlots;
}


#endif
