#ifndef AliReaderESD_H
#define AliReaderESD_H
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

#include "AliReader.h"
#include <TString.h>
class TFile;
class AliRunLoader;
class AliESD;
class AliESDtrack;

class AliReaderESD: public AliReader
{
  public:
    AliReaderESD(const Char_t* esdfilename = "AliESDs.root", const Char_t* galfilename = "galice.root");

    AliReaderESD(TObjArray* dirs,const Char_t* esdfilename = "AliESDs.root", const Char_t* galfilename = "galice.root");

    virtual ~AliReaderESD();
    
    void          Rewind();
    
    void          ReadSimulatedData(Bool_t flag){fReadSim = flag;}//switches reading MC data
    Bool_t        ReadsRec() const {return kTRUE;}
    Bool_t        ReadsSim() const {return fReadSim;}
    void          SetCheckParticlePID(Bool_t flag){fCheckParticlePID = flag;}
    void          SetReadMostProbableOnly(Bool_t flag){fReadMostProbableOnly = flag;}
            
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
    void          SetITSTrackPoints(Bool_t flag = kTRUE){fITSTrackPoints = flag;}
    void          MustTPC(Bool_t flag){fMustTPC = flag;}

    void          SetReadCentralBarrel(Bool_t flag){fReadCentralBarrel = flag;}
    void          SetReadMuon(Bool_t flag){fReadMuon = flag;}
    void          SetReadPHOS(Bool_t flag){fReadPHOS = flag;}

    enum ESpecies {kESDElectron = 0, kESDMuon, kESDPion, kESDKaon, kESDProton, kNSpecies};
    static Int_t  GetSpeciesPdgCode(ESpecies spec);
    
    Int_t         ReadESD(AliESD* esd);
    Int_t         ReadESDCentral(AliESD* esd);
    Int_t         ReadESDMuon(AliESD* esd);
    Int_t         ReadESDPHOS(AliESD* /*esd*/){return 0;}

  protected:
    virtual Int_t         ReadNext();
 
    virtual TFile*        OpenFile(Int_t evno);//opens files to be read for given event

    Bool_t        CheckTrack(AliESDtrack* t) const;
    
    TString       fESDFileName;//name of the file with tracks
    TString       fGAlFileName;//name of the file with tracks
    TFile*        fFile;//! pointer to current ESD file
    AliRunLoader* fRunLoader;//!Run Loader
    TIter*        fKeyIterator;//!iterator over keys in ESD file
    Bool_t        fReadSim;//flag indicating wether to read particles from kinematics
    Bool_t        fCheckParticlePID;//flag indicating to perform the check on PID of simulated particle - usefull in resoluion analysis
    Bool_t        fReadMostProbableOnly;//flag indicating to read ony one incarnation with the highest probability
    Int_t         fNTrackPoints;//number of track points; if==0 track points are not created
    Float_t       fdR;//spacing between points (along radius) in cm
                      //Track Points are needed for Anti-Merging Cut
    
    Bool_t        fClusterMap;//Flag indicating if Claster Map should be created for each track
                              //Claster map is needed for Anti-Splitting Cut

    Bool_t        fITSTrackPoints;//Flag indicalting if track positions in ITS are to be read
                                  //currently we use only position at first pixels wich are
	              //used by anti-merging cut in non-id analysis

    Bool_t        fMustTPC;// must be reconstructed in TPC -> reject tracks reconstructed ITS stand alone

    Bool_t        fReadCentralBarrel; // Flag for reading ESD central track 
    Bool_t        fReadMuon;// Flag for reading ESD Muon track 
    Bool_t        fReadPHOS;// Flag for reading ESD Phos 

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
    ClassDef(AliReaderESD,1)
};


#endif
