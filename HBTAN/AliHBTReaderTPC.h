#ifndef AliHBTReaderTPC_H
#define AliHBTReaderTPC_H
//______________________________________________
//
// class AliHBTReaderTPC
//
// reader for TPC tracks
// needs galice.root
// 
// more info: http://aliweb.cern.ch/people/skowron/analyzer/index.html
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////////////////////

#include "AliHBTReader.h"
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
    AliHBTReaderTPC(const AliHBTReaderTPC& in);
    
    virtual ~AliHBTReaderTPC();
    
    AliHBTReaderTPC& operator=(const AliHBTReaderTPC& in);
    
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
    void          SetClusterMap(Bool_t flag = kTRUE){fClusterMap = flag;}

  protected:
    Int_t         ReadNext();
    Int_t         OpenNextSession();
    void          DoOpenError(const char* msgfmt, ...);
    
    TString       fFileName;//name of the file with galice.root
    AliRunLoader* fRunLoader;//!RL
    AliTPCLoader* fTPCLoader;//!TPCLoader
    Float_t       fMagneticField;//magnetic field value that was enforced while reading
    Bool_t        fUseMagFFromRun;//flag indicating if using field specified in gAlice (kTRUE)
                               // or enforece other defined by fMagneticField
    
    Int_t         fNTrackPoints;//number of track points; if==0 track points are not created
    Float_t       fdR;//spacing between points (along radius) in cm
                      //Track Points are needed for Anti-Merging Cut
    
    Bool_t        fClusterMap;//Flag indicating if Claster Map should be created for each track
                              //Claster map is needed for Anti-Splitting Cut

    //Cut Parameters specific to TPC tracks
        
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
    
    Bool_t CheckTrack(AliTPCtrack* t) const;

    ClassDef(AliHBTReaderTPC,3)
};


#endif
