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


class AliLHCData : public TObject
{
 public:
  enum          {kStart,kNStor};
  enum BeamID_t {kBeam1,kBeam2};
  enum Proj_t   {kX,kY};
  enum Side_t   {kLeft,kRight};
  enum Collim_t {kTCTVB4L2, kTCTVB4R2, kTCLIA4R2, kNCollimators};
  enum ColJaw_t {kGapDn,kGapUp,kLeftDn,kLeftUp,kRightDn,kRightUp,kNJaws};
  enum          {kMaxBSlots = 3564};
  enum          {kMarginSOR = 60*60*24*30, // use margin of 30 days for SOR, when looking for the 1st record
		 kMarginEOR = 60*15};      // use margin of 15 min for EOR, when looking for the last record
  //
  enum {kIntTot,kIntTotAv,kIntBunchAv,
	kLumAcqMode,kLumTot,kLumTotErr,kLumBunch,kLumBunchErr,kLumCrossAng,kLumCrossAngErr,
	kBunchConf,kFillNum,kBunchLgtNB,kBunchLgt,kBunchLgtFillB,
	kRCInjSch,kRCBeta,kRCCrossAng,kRCVang,
	kBeamSzAcqMode,kBeamSzSigH,kBeamSzSigV,kBeamSzEmittH,kBeamSzEmittV,kBeamSzSigHErr,kBeamSzSigVErr,
	kCollPos};
  //
  //le
 public:
  //
 AliLHCData() : fTMin(0),fTMax(1e10),fFillNumber(0),fData(0),fFile2Process(0),fMap2Process(0) {Clear();}
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
  //
  virtual void          Print(const Option_t *opt="")                const;
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
  void           FlagInteractingBunches(const Int_t beam1[2],const Int_t beam2[2]);
  TObject*       FindRecValidFor(int start,int nrec, double tstamp) const;
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
  const Int_t* GetOffsCollimatorJawPos(int coll,int jaw)  const;
  //
  const TObjArray&  GetData()                             const {return fData;}
  //
  // analysis methods
  Int_t GetMeanIntensity(int beamID, Double_t &colliding, Double_t &noncolliding) const;
  //
 protected:
  //
  Bool_t                FillData(double tmin=0, double tmax=1.e20);
  virtual void          Clear(const Option_t *opt="");
  void                  PrintAux(Bool_t full,const Int_t refs[2],const Option_t *opt="") const;
  TObjArray*            GetDCSEntry(const char* key,int &entry,int &last,double tmin,double tmax) const;
  Int_t                 FillScalarRecord(  int refs[2], const char* rec, const char* recErr=0);
  Int_t                 FillBunchConfig(   int refs[2], const char* rec);
  Int_t                 FillStringRecord(  int refs[2], const char* rec);
  Int_t                 FillAcqMode(       int refs[2], const char* rec);
  Int_t                 FillBunchInfo(     int refs[2], const char* rec,int ibm, Bool_t inRealSlots);
  Int_t                 FillBCLuminosities(int refs[2], const char* rec, const char* recErr, Int_t useBeam);
  //
  Int_t                 ExtractInt(AliDCSArray* dcsArray,Int_t el)    const;
  Double_t              ExtractDouble(AliDCSArray* dcsArray,Int_t el) const;
  TString&              ExtractString(AliDCSArray* dcsArray)          const;
 AliLHCData(const AliLHCData& src) : TObject(src),fTMin(0),fTMax(0),fFillNumber(0),fData(0),fFile2Process(0),fMap2Process(0) { /*dummy*/ }
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
  Int_t           fLuminAcqMode[2][2];                // luminosity acquisition mode                   : Int    | kLumAcqMode
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
  const Char_t*   fFile2Process;                      //! name of DCS file
  const TMap*     fMap2Process;                       //! DCS map to process 

  ClassDef(AliLHCData,1)
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


#endif
