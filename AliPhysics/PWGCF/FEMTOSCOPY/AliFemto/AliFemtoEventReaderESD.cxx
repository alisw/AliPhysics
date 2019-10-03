///
/// \file AliFemto/AliFemtoEventReaderESD.cxx
///

#include "AliFemtoEventReaderESD.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"

//#include "TSystem.h"

#include "AliFmPhysicalHelixD.h"
#include "AliFmThreeVectorF.h"

#include "SystemOfUnits.h"

#include "AliFemtoEvent.h"
#include "AliFemtoModelHiddenInfo.h"

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoEventReaderESD);
  /// \endcond
#endif


#if !(ST_NO_NAMESPACES)
  using namespace units;
#endif

using namespace std;
//____________________________
//constructor with 0 parameters , look at default settings
AliFemtoEventReaderESD::AliFemtoEventReaderESD():
  fInputFile(" "),
  fFileName(" "),
  fConstrained(true),
  fReadInner(false),
  fNumberofEvent(0),
  fCurEvent(0),
  fTree(nullptr),
  fEsdFile(nullptr),
  fEvent(nullptr)
{
  // default constructor
}

AliFemtoEventReaderESD::AliFemtoEventReaderESD(const AliFemtoEventReaderESD &aReader) :
  AliFemtoEventReader(aReader),
  fInputFile(aReader.fInputFile),
  fFileName(aReader.fFileName),
  fConstrained(aReader.fConstrained),
  fReadInner(aReader.fReadInner),
  fNumberofEvent(aReader.fNumberofEvent),
  fCurEvent(aReader.fCurEvent),
  fTree(nullptr),
  fEsdFile(new TFile(aReader.fEsdFile->GetName())),
  fEvent(new AliESDEvent())
{
  // copy constructor
  //  fTree = aReader.fTree->CloneTree();
  //  fEvent = new AliESD(*aReader.fEvent);
}
//__________________
//Destructor
AliFemtoEventReaderESD::~AliFemtoEventReaderESD()
{
  // destructor
  //delete fListOfFiles;
  delete fTree;
  delete fEvent;
  delete fEsdFile;
}

//__________________
AliFemtoEventReaderESD& AliFemtoEventReaderESD::operator=(const AliFemtoEventReaderESD &aReader)
{
  // assignment operator
  if (this == &aReader)
    return *this;

  AliFemtoEventReader::operator=(aReader);

  fInputFile = aReader.fInputFile;
  fFileName  = aReader.fFileName;
  fConstrained = aReader.fConstrained;
  fReadInner = aReader.fReadInner;
  fNumberofEvent = aReader.fNumberofEvent;
  fCurEvent = aReader.fCurEvent;

  delete fTree;
  fTree = nullptr;
  //  fTree = aReader.fTree->CloneTree();

  delete fEvent;
  fEvent = new AliESDEvent();

  delete fEsdFile;
  fEsdFile = new TFile(aReader.fEsdFile->GetName());

  return *this;
}
//__________________
AliFemtoString AliFemtoEventReaderESD::Report()
{
  // create reader report
  AliFemtoString temp = "\n This is the AliFemtoEventReaderESD\n";
  return temp;
}

//__________________
void AliFemtoEventReaderESD::SetInputFile(const char* inputFile)
{
  //setting the name of file where names of ESD file are written
  //it takes only this files which have good trees
  fInputFile = inputFile;
  cout<<"Input File set on "<<fInputFile<<endl;

  delete fTree;
  fTree = new TChain("esdTree");

  ifstream infile(inputFile);

  if (infile.good()) {
    // adding all given files with good tree inside
    std::string buffer;
    while (infile.eof() == false) {
      std::getline(infile, buffer);
      TFile *esdFile=TFile::Open(buffer.c_str(), "READ");
      if (esdFile) {
        TTree* tree = static_cast<TTree*>(esdFile->Get("esdTree"));
        if (tree) {
          cout<<"putting file  "<<buffer<<" into analysis"<<endl;
          fTree->AddFile(esdFile->GetName());
          delete tree;
        }
        esdFile->Close();
      }
      delete esdFile;
    }
  }
}

void AliFemtoEventReaderESD::SetConstrained(const bool constrained)
{
  fConstrained=constrained;
}

bool AliFemtoEventReaderESD::GetConstrained() const
{
  return fConstrained;
}

void AliFemtoEventReaderESD::SetReadTPCInner(const bool readinner)
{
  fReadInner=readinner;
}

bool AliFemtoEventReaderESD::GetReadTPCInner() const
{
  return fReadInner;
}

AliFemtoEvent* AliFemtoEventReaderESD::ReturnHbtEvent()
{
  // read in a next hbt event from the chain
  // convert it to AliFemtoEvent and return
  // for further analysis
  AliFemtoEvent *hbtEvent = 0;

  if (fCurEvent==fNumberofEvent) { //open next file
    if (fNumberofEvent==0) {
      //	  delete fEvent;//added 1.04.2007
      fEvent=new AliESDEvent();
      //	  delete fTree;
      //fTree=0;
      //	  delete fEsdFile;

      //ESD data
      //	  fEsdFile=TFile::Open(fFileName.c_str(),"READ");
      //	  fTree = (TTree*) fEsdFile->Get("esdTree");
      //	  fTree->SetBranchAddress("ESD", &fEvent);
      fTree->SetBranchStatus("MuonTracks*",0);
      fTree->SetBranchStatus("PmdTracks*",0);
      fTree->SetBranchStatus("TrdTracks*",0);
      fTree->SetBranchStatus("V0s*",0);
      fTree->SetBranchStatus("Cascades*",0);
      fTree->SetBranchStatus("Kinks*",0);
      fTree->SetBranchStatus("CaloClusters*",0);
      fTree->SetBranchStatus("AliRawDataErrorLogs*",0);
      fTree->SetBranchStatus("ESDfriend*",0);
      fEvent->ReadFromTree(fTree);

      fNumberofEvent=fTree->GetEntries();
      cout<<"Number of Entries in file "<<fNumberofEvent<<endl;
      fCurEvent=0;
      //sim data
    } else { //no more data to read
      cout<<"no more files "<<hbtEvent<<endl;
      fReaderStatus=1;
      return hbtEvent;
    }
  }

  cout<<"starting to read event "<<fCurEvent<<endl;
  fTree->GetEvent(fCurEvent);//getting next event
  cout << "Read event " << fEvent << " from file " << fTree << endl;
  //  vector<int> tLabelTable;//to check labels

  hbtEvent = new AliFemtoEvent;
  //setting basic things
  //  hbtEvent->SetEventNumber(fEvent->GetEventNumber());
  hbtEvent->SetRunNumber(fEvent->GetRunNumber());
  //hbtEvent->SetNumberOfTracks(fEvent->GetNumberOfTracks());
  hbtEvent->SetMagneticField(fEvent->GetMagneticField()*kilogauss);//to check if here is ok
  hbtEvent->SetZDCN1Energy(fEvent->GetZDCN1Energy());
  hbtEvent->SetZDCP1Energy(fEvent->GetZDCP1Energy());
  hbtEvent->SetZDCN2Energy(fEvent->GetZDCN2Energy());
  hbtEvent->SetZDCP2Energy(fEvent->GetZDCP2Energy());
  hbtEvent->SetZDCEMEnergy(fEvent->GetZDCEMEnergy());
  hbtEvent->SetZDCParticipants(fEvent->GetZDCParticipants());
  hbtEvent->SetTriggerMask(fEvent->GetTriggerMask());
  hbtEvent->SetTriggerCluster(fEvent->GetTriggerCluster());

  //Vertex
  double fV1[3];
  fEvent->GetVertex()->GetXYZ(fV1);

  AliFmThreeVectorF vertex(fV1[0],fV1[1],fV1[2]);
  hbtEvent->SetPrimVertPos(vertex);

  //starting to reading tracks
  int nofTracks=0;  //number of reconstructed tracks in event
  nofTracks=fEvent->GetNumberOfTracks();
  int realnofTracks=0;//number of track which we use ina analysis
  cout << "Event has " << nofTracks << " tracks " << endl;

  for (int i=0;i<nofTracks;i++) {
    bool  tGoodMomentum=true; //flaga to chcek if we can read momentum of this track

    AliFemtoTrack* trackCopy = new AliFemtoTrack();
    const AliESDtrack *esdtrack=fEvent->GetTrack(i);//getting next track
    //      const AliESDfriendTrack *tESDfriendTrack = esdtrack->GetFriendTrack();

    trackCopy->SetCharge((short)esdtrack->GetSign());

    //in aliroot we have AliPID
    //0-electron 1-muon 2-pion 3-kaon 4-proton 5-photon 6-pi0 7-neutron 8-kaon0 9-eleCon
    //we use only 5 first
    double esdpid[5];
    esdtrack->GetESDpid(esdpid);
    trackCopy->SetPidProbElectron(esdpid[0]);
    trackCopy->SetPidProbMuon(esdpid[1]);
    trackCopy->SetPidProbPion(esdpid[2]);
    trackCopy->SetPidProbKaon(esdpid[3]);
    trackCopy->SetPidProbProton(esdpid[4]);

    double pxyz[3];
    if (fReadInner == true) {
      if (auto *inner_param = esdtrack->GetTPCInnerParam()) {
        AliExternalTrackParam param(*inner_param);
        param.PropagateToDCA(fEvent->GetPrimaryVertex(), fEvent->GetMagneticField(), 10000);
        param.GetPxPyPz(pxyz);//reading noconstarined momentum

        AliFemtoModelHiddenInfo *tInfo = new AliFemtoModelHiddenInfo();
        tInfo->SetPDGPid(211);
        tInfo->SetTrueMomentum(pxyz[0], pxyz[1], pxyz[2]);
        tInfo->SetMass(0.13957);
        trackCopy->SetHiddenInfo(tInfo);
      }
    }

    tGoodMomentum = fConstrained
                  ? esdtrack->GetConstrainedPxPyPz(pxyz)  // reading constrained momentum
                  : esdtrack->GetPxPyPz(pxyz);  //reading noconstarined momentum

    AliFemtoThreeVector v(pxyz[0],pxyz[1],pxyz[2]);
    trackCopy->SetP(v);//setting momentum
    trackCopy->SetPt(sqrt(pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]));
    const AliFmThreeVectorD ktP(pxyz[0],pxyz[1],pxyz[2]);
    if (ktP.Mag() == 0) {
      delete trackCopy;
      continue;
    }

    const AliFmThreeVectorD origin(fV1[0],fV1[1],fV1[2]);
    //setting helix I do not if it is ok
    AliFmPhysicalHelixD helix(ktP,origin,(double)(fEvent->GetMagneticField())*kilogauss,(double)(trackCopy->Charge()));
    trackCopy->SetHelix(helix);

    trackCopy->SetTrackId(esdtrack->GetID());
    trackCopy->SetFlags(esdtrack->GetStatus());
    trackCopy->SetLabel(esdtrack->GetLabel());

    //some stuff which could be useful
    float impact[2];
    float covimpact[3];
    esdtrack->GetImpactParameters(impact,covimpact);
    trackCopy->SetImpactD(impact[0]);
    trackCopy->SetImpactZ(impact[1]);
    trackCopy->SetCdd(covimpact[0]);
    trackCopy->SetCdz(covimpact[1]);
    trackCopy->SetCzz(covimpact[2]);
    trackCopy->SetITSchi2(esdtrack->GetITSchi2());
    trackCopy->SetITSncls(esdtrack->GetNcls(0));
    trackCopy->SetTPCchi2(esdtrack->GetTPCchi2());
    trackCopy->SetTPCncls(esdtrack->GetTPCNcls());
    trackCopy->SetTPCnclsF(esdtrack->GetTPCNclsF());
    trackCopy->SetTPCsignalN((short)esdtrack->GetTPCsignalN()); //due to bug in aliesdtrack class
    trackCopy->SetTPCsignalS(esdtrack->GetTPCsignalSigma());

    trackCopy->SetTPCClusterMap(esdtrack->GetTPCClusterMap());
    trackCopy->SetTPCSharedMap(esdtrack->GetTPCSharedMap());

    double pvrt[3];
    fEvent->GetPrimaryVertex()->GetXYZ(pvrt);

    double xtpc[3];
    esdtrack->GetInnerXYZ(xtpc);
    xtpc[2] -= pvrt[2];
    trackCopy->SetNominalTPCEntrancePoint(xtpc);

    esdtrack->GetOuterXYZ(xtpc);
    xtpc[2] -= pvrt[2];
    trackCopy->SetNominalTPCExitPoint(xtpc);

    int indexes[3];
    for (int ik=0; ik<3; ik++) {
      indexes[ik] = esdtrack->GetKinkIndex(ik);
    }

    trackCopy->SetKinkIndexes(indexes);
    //decision if we want this track
    //if we using diffrent labels we want that this label was use for first time
    //if we use hidden info we want to have match between sim data and ESD
    if (tGoodMomentum==true) {
      hbtEvent->TrackCollection()->push_back(trackCopy);//adding track to analysis
      realnofTracks++;//real number of tracks
    } else {
      delete  trackCopy;
    }
  }

  hbtEvent->SetNumberOfTracks(realnofTracks);//setting number of track which we read in event
  fCurEvent++;
  cout<<"end of reading nt "<<nofTracks<<" real number "<<realnofTracks<<endl;
//   if (fCurEvent== fNumberofEvent)//if end of current file close all
//     {
//       fTree->Reset();
//       delete fTree;
//       fEsdFile->Close();
//     }
  return hbtEvent;
}
