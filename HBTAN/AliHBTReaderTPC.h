#ifndef AliHBTReaderTPC_H
#define AliHBTReaderTPC_H

#include "AliHBTReader.h"

//Multi file reader for TPC
//
//This reader reads tracks AliTPCtracks.root
//                  particles form gAlice
//Piotr.Skowronski@cern.ch
//more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html

#include <TString.h>
class TFile;
class TArrayF;
class AliRunLoader;
class AliTPCLoader;
class AliTPCtrack;

class AliHBTReaderTPC: public AliHBTReader
{
  public:
    AliHBTReaderTPC();
    AliHBTReaderTPC(const Char_t* galicefilename);
    AliHBTReaderTPC(TObjArray* dirs, const Char_t* galicefilename = "galice.root");

    virtual ~AliHBTReaderTPC();
    
    void          Rewind();
    
    Bool_t        ReadsTracks() const {return kTRUE;}
    Bool_t        ReadsParticles() const {return kTRUE;}
    
    void          SetMagneticField(Float_t mf){fMagneticField=mf;}
    void          UseMagneticFieldFromRun(Bool_t flag = kTRUE){fUseMagFFromRun=flag;}
    
    void          SetNClustersRange(Int_t min,Int_t max);
    void          SetChi2PerCluserRange(Float_t min, Float_t max);
    void          SetC00Range(Float_t min, Float_t max);
    void          SetC11Range(Float_t min, Float_t max);
    void          SetC22Range(Float_t min, Float_t max);
    void          SetC33Range(Float_t min, Float_t max);
    void          SetC44Range(Float_t min, Float_t max);
    void          SetNumberOfTrackPoints(Int_t n = 5,Float_t dr = 30.0) {fNTrackPoints = n; fdR = dr;}
    Int_t         GetNumberOfTrackPoints() const {return fNTrackPoints;}
  protected:
    //in the future this class is will read global tracking
    Int_t         ReadNext();
    Int_t         OpenNextSession();
    void          DoOpenError(const char* msgfmt, ...);
    
    TString       fFileName;//name of the file with galice.root
    AliRunLoader* fRunLoader;//!RL
    AliTPCLoader* fTPCLoader;//!TPCLoader
    Float_t       fMagneticField;//magnetic field value that was enforced while reading
    Bool_t        fUseMagFFromRun;//flag indicating if using field specified in gAlice (kTRUE)
                               // or enforece other defined by fMagneticField
    
    Int_t         fNTrackPoints;//number of track points
    Float_t       fdR;//spacing between points (along radius) in cm
        
    Int_t         fNClustMin;//Number of clusters min value
    Int_t         fNClustMax;//Number of clusters max value
    
    Float_t       fNChi2PerClustMin;//Chi^2 per number of clusters min value
    Float_t       fNChi2PerClustMax;//Chi^2 per number of clusters max value

    Float_t       fC00Min;//C00 (0th diagonal element of covariance matrix) min value
    Float_t       fC00Max;//C00 (0th diagonal element of covariance matrix) max value
            
    Float_t       fC11Min;//C11 (1th diagonal element of covariance matrix) min value
    Float_t       fC11Max;//C11 (1th diagonal element of covariance matrix) max value
    
    Float_t       fC22Min;//C22 (2th diagonal element of covariance matrix) min value
    Float_t       fC22Max;//C22 (2th diagonal element of covariance matrix) max value
    
    Float_t       fC33Min;//C33 (3th diagonal element of covariance matrix) min value
    Float_t       fC33Max;//C33 (3th diagonal element of covariance matrix) max value
    
    Float_t       fC44Min;//C44 (4th diagonal element of covariance matrix) min value
    Float_t       fC44Max;//C44 (4th diagonal element of covariance matrix) max value

  private:
    
    Bool_t CheckTrack(AliTPCtrack* t);
  public:
    ClassDef(AliHBTReaderTPC,3)
};


#endif
