#ifndef ALIHBTREADERESD_H
#define ALIHBTREADERESD_H
//___________________________________________________________________________
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Multi file reader for ESD                                               //
//                                                                         //
// This reader reads tracks from Event Summary Data                        //
// do not read particles                                                   //
// Piotr.Skowronski@cern.ch                                                //
// more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html    //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "AliHBTReader.h"
#include <TString.h>
class TFile;
class AliRunLoader;
class AliESD;
class AliESDtrack;

class AliHBTReaderESD: public AliHBTReader
{
  public:
    AliHBTReaderESD(const Char_t* esdfilename = "AliESDs.root", const Char_t* galfilename = "galice.root");

    AliHBTReaderESD(TObjArray* dirs,const Char_t* esdfilename = "AliESDs.root", const Char_t* galfilename = "galice.root");

    virtual ~AliHBTReaderESD();
    
    void          Rewind();
    
    void          ReadParticles(Bool_t flag){fReadParticles = flag;}
    Bool_t        ReadsTracks() const {return kTRUE;}
    Bool_t        ReadsParticles() const {return fReadParticles;}
    void          SetCheckParticlePID(Bool_t flag){fCheckParticlePID = flag;}
        
    void          ReadDataTPC(){}
    void          ReadDataITS(){}

    void          SetTPCNClustersRange(Int_t min,Int_t max);
    void          SetTPCChi2PerCluserRange(Float_t min, Float_t max);
    
    void          SetChi2Range(Float_t min, Float_t max);
    void          SetC00Range(Float_t min, Float_t max);
    void          SetC11Range(Float_t min, Float_t max);
    void          SetC22Range(Float_t min, Float_t max);
    void          SetC33Range(Float_t min, Float_t max);
    void          SetC44Range(Float_t min, Float_t max);
    void          SetNumberOfTrackPoints(Int_t n = 5,Float_t dr = 30.0) {fNTrackPoints = n; fdR = dr;}
    Int_t         GetNumberOfTrackPoints() const {return fNTrackPoints;}
    void          SetClusterMap(Bool_t flag = kTRUE){fClusterMap = flag;}
    void          MustTPC(Bool_t flag){fMustTPC = flag;}
    
    enum ESpecies {kESDElectron = 0, kESDMuon, kESDPion, kESDKaon, kESDProton, kNSpecies};
    static Int_t  GetSpeciesPdgCode(ESpecies spec);//skowron
    
    Int_t         ReadESD(AliESD* esd);
    
  protected:
    Int_t         ReadNext();
    TFile*        OpenFile(Int_t evno);//opens files to be read for given event
    Bool_t        CheckTrack(AliESDtrack* t) const;
    
    TString       fESDFileName;//name of the file with tracks
    TString       fGAlFileName;//name of the file with tracks
    TFile*        fFile;//! pointer to current ESD file
    AliRunLoader* fRunLoader;//!Run Loader
    TIter*        fKeyIterator;
    Bool_t        fReadParticles;//flag indicating wether to read particles from kinematics
    Bool_t        fCheckParticlePID;//flag indicating to perform the check on PID of simulated particle
    
    Int_t         fNTrackPoints;//number of track points; if==0 track points are not created
    Float_t       fdR;//spacing between points (along radius) in cm
                      //Track Points are needed for Anti-Merging Cut
    
    Bool_t        fClusterMap;//Flag indicating if Claster Map should be created for each track
                              //Claster map is needed for Anti-Splitting Cut

    Bool_t        fMustTPC;// must be reconstructed in TPC -> reject tracks reconstructed ITS stand alone
    
    //Cut Parameters specific to TPC tracks
        
    Int_t         fNTPCClustMin;//Number of clusters min value
    Int_t         fNTPCClustMax;//Number of clusters max value
    
    Float_t       fTPCChi2PerClustMin;//Chi^2 per number of clusters min value
    Float_t       fTPCChi2PerClustMax;//Chi^2 per number of clusters max value


    // Required parameters at vertex
    Float_t       fChi2Min;//Chi^2 min value
    Float_t       fChi2Max;//Chi^2 max value

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

    // Required parameters at TPC Inner Layer
    Float_t       fTPCC00Min;//C00 (0th diagonal element of covariance matrix) min value
    Float_t       fTPCC00Max;//C00 (0th diagonal element of covariance matrix) max value
            
    Float_t       fTPCC11Min;//C11 (1th diagonal element of covariance matrix) min value
    Float_t       fTPCC11Max;//C11 (1th diagonal element of covariance matrix) max value
    
    Float_t       fTPCC22Min;//C22 (2th diagonal element of covariance matrix) min value
    Float_t       fTPCC22Max;//C22 (2th diagonal element of covariance matrix) max value
    
    Float_t       fTPCC33Min;//C33 (3th diagonal element of covariance matrix) min value
    Float_t       fTPCC33Max;//C33 (3th diagonal element of covariance matrix) max value
    
    Float_t       fTPCC44Min;//C44 (4th diagonal element of covariance matrix) min value
    Float_t       fTPCC44Max;//C44 (4th diagonal element of covariance matrix) max value
    
  private:
    ClassDef(AliHBTReaderESD,2)
};


#endif
