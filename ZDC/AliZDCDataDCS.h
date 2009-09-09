#ifndef AliZDCDataDCS_H
#define AliZDCDataDCS_H

////////////////////////////////////////////////
//  Class for ZDC DCS data                    //
////////////////////////////////////////////////


#include <TMap.h>

class AliZDCDataDCS : public TObject {
public:
    enum {kNAliases=28};

    AliZDCDataDCS();
    AliZDCDataDCS(Int_t nRun, UInt_t startTime, UInt_t endTime, 
                  UInt_t startTimeDCSQuery, UInt_t endTimeDCSQuery);
    AliZDCDataDCS(const AliZDCDataDCS & data);
    AliZDCDataDCS& operator=(const AliZDCDataDCS & data);
    ~AliZDCDataDCS();

    void SetRun(Int_t run) {fRun = run;}
    void SetStartTime(Int_t startTime) {fStartTime = startTime;}
    void SetEndTime(Int_t endTime) {fEndTime = endTime;}
    void SetStartTimeDCSQuery(Int_t startTimeDCSQuery) {fStartTimeDCSQuery = startTimeDCSQuery;}
    void SetEndTimeDCSQuery(Int_t endTimeDCSQuery) {fEndTimeDCSQuery = endTimeDCSQuery;}
    void SetCalibData(Float_t *val) {for(Int_t i=0; i<kNAliases; i++) fCalibData[i] = val[i];}
    void SetCalibData(Int_t i, Float_t val) {fCalibData[i] = val;} 
    //
    Int_t GetRun() const {return fRun;}
    Int_t GetStartTime() const {return fStartTime;}
    Int_t GetEndTime() const {return fEndTime;}
    Int_t GetStartTimeDCSQuery() const {return fStartTimeDCSQuery;}
    Int_t GetEndTimeDCSQuery() const {return fEndTimeDCSQuery;}
    Float_t GetCalibData() const {return *fCalibData;}
    Float_t GetCalibData(Int_t i) const {return fCalibData[i];}
    Float_t GetTimeStamp() const {return *fTimeStamp;}
    Float_t GetTimeStamp(Int_t i) const {return fTimeStamp[i];}

    Bool_t ProcessData(TMap& aliasMap);

    const char* GetAliasName(UInt_t pos) const
    		{return pos<kNAliases ? fAliasNames[pos].Data() : 0;}

private:
    void Init();
    void Introduce(UInt_t numAlias, const TObjArray* aliasArr);

    Int_t  fRun; 	    // Run number
    UInt_t fStartTime;      // Start of run time
    UInt_t fEndTime;	    // End of run time
    Int_t  fStartTimeDCSQuery; // start time DCSQuery
    Int_t  fEndTimeDCSQuery;   // end time DCSQuery

    TString  fAliasNames[kNAliases];// Name of the aliases provided by the DCS
    Double_t fCalibData[kNAliases]; // Array containing calibration data
    Float_t  fTimeStamp[kNAliases]; // Array containing time stamps

    Bool_t fIsProcessed;	    // Flag set when data are processed

    ClassDef(AliZDCDataDCS, 5);
};

#endif
