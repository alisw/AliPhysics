#ifndef ALIREADERAOD_H
#define ALIREADERAOD_H
//______________________________________________________________________________
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class AliReaderAOD                                                         //
//                                                                            //
// Reader and Writer for AOD format.                                          //
// AODs are stored in a tree named by the variable fgkTreeName.               //
// There is stored 1 or 2 branches. Each of them stores AOD objects           //
// First branch is named by the variable fgkReconstructedDataBranchName       //
// ("reconstructed.") and keeps reconstructed data.                           //
// Second branch is called by the variable fgkSimulatedDataBranchName         //
// ("simulated.") and stores Monte carlo truth. If both branches are present  //
// AODs are parallel, i.e. nth particle in one branch corresponds to the nth  //
// particle in the other one.                                                 //
//                                                                            //
// Since we accept different formats of particles that are stored in AODs     //
// reader must take care of that fact: clean buffer if the next file contains //
// different particle type.                                                   //
//                                                                            //
// If no cuts are specified in a reader, it reuturns pointer to the           //
// buffers. In the other case data are copied to the onother AOD (filtering   //
// out particles that do not pass a cut), thus reading is slower.             //
//                                                                            //
// Piotr.Skowronski@cern.ch                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////



#include "AliReader.h"

class TTree;
class TFile;

class AliReaderAOD: public AliReader
{
  public:
    AliReaderAOD(const Char_t* aodfilename = "AOD.root");
    virtual ~AliReaderAOD();

    void          ReadSimulatedData(Bool_t flag){fReadSim = flag;}//switches reading MC data
    void          ReadReconsructedData(Bool_t flag){fReadRec = flag;}//switches reading MC data
    Bool_t        ReadsRec() const {return fReadRec;}
    Bool_t        ReadsSim() const {return fReadSim;}

    void          Rewind();


    static Int_t WriteAOD(AliReader* reader, const char* outfilename = "AliAOD.root", //reads tracks from runs and writes them to file
                          const char* pclassname = "AliAODParticle", Bool_t multcheck = kFALSE);
    
  protected:
    virtual Int_t         ReadNext();
    virtual Int_t         OpenFile(Int_t evno);//opens files to be read for given event
    
    virtual Int_t         ReadRecAndSim();
    virtual Int_t         ReadRec();
    virtual Int_t         ReadSim();
    
    static const TString  fgkTreeName;//name of branch holding simulated data
    static const TString  fgkReconstructedDataBranchName;//name of branch holding reconstructed data
    static const TString  fgkSimulatedDataBranchName;//name of branch holding simulated data
    
  
  private:
    TString fFileName;//File name
    
    Bool_t  fReadSim;//indicates if to read simulated data
    Bool_t  fReadRec;//indicates if to read simulated data

    TTree*        fTree;//!tree
    TFile*        fFile;//!file
    AliAOD*       fSimBuffer;//!buffer array that tree is read to
    AliAOD*       fRecBuffer;//!
    
    
    ClassDef(AliReaderAOD,1)
};

#endif
