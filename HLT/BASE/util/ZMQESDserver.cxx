#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliHLTZMQhelpers.h"
#include "AliOptionParser.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "signal.h"
#include "zmq.h"
#include <iostream>
#include <memory>
#include <pthread.h>
#include <string>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

using namespace AliZMQhelpers;

// methods
Int_t ProcessOptionString(TString arguments, Bool_t verbose = kFALSE);
Int_t InitZMQ();
Int_t Run();

Int_t HandleRequest(aliZMQmsg::iterator block, void* /*socket*/ = NULL);
Int_t DoSend(void* socket);
Int_t HandleControlSequence(aliZMQmsg::iterator block, void* socket);

void SetLastPushBackTime();

template <typename T> class o2vec;
template <typename T>
void alizmq_deleteArray(void* data, void*);
template <typename T>
int alizmq_msg_add_array(aliZMQmsg* message, o2vec<T>&& vec);

// configuration vars
Bool_t fVerbose = kFALSE;
TString fZMQconfigOUT = "PUSH";
Int_t fZMQmaxQueueSize = 10;
Int_t fZMQtimeout = -1;
int fZMQlingerTime = 30000; //30s default linger

std::vector<std::string> fFilesIn;

std::string fInfo; // cache for the info string

std::string fID; // merger ID/name

long fPushbackPeriod = -1;    // in seconds, -1 means never
long fRequestPeriod = -1;     // in seconds, -1 means never
long fRequestTimeout = 10000; // default 10s
aliZMQrootStreamerInfo fSchema;
bool fSchemaOnRequest = false;
bool fSchemaOnSend = false;
int fCompression = 1;
bool fgTerminationSignaled = false;
bool fOnce = true;

unsigned long long fLastPushBackTime = 0;
unsigned long long fLastRequestTime = 0;
struct timeval fCurrentTimeStruct;

uint64_t fBytesOUT = 0;
unsigned long fSecondsStart = 0;
unsigned long fSecondsStop = 0;
unsigned long long fNumberOfMessagesSent = 0;

// ZMQ stuff
void* fZMQcontext = NULL; // ze zmq context

void* fZMQout = NULL; // the monitoring socket, here we publish a copy of the data
void* fZMQresetBroadcast = NULL;
bool fReconfigureZMQ = true;

const char* fUSAGE = "ZMQESDserver options: \n"
                     " -id : some string identifier\n"
                     " -in : comma separated list of input files, e.g. AliESDs.root,../AliESDs.root\n"
                     " -out : data out, format is zmq config string, e.g. PUSH+tcp://localhost:123123 or REP@tcp://*:123123\n"
                     " -once : [0|1] run once, or run forever"
                     " -ZMQlinger : [ms] - how long to wait for socket to send pending data before forcing close"
                     " -Verbose : print some info\n";
//_______________________________________________________________________________________
void sig_handler(int signo)
{
  if (signo == SIGINT)
    printf(" received SIGINT\n");
  fgTerminationSignaled = true;
}

//_______________________________________________________________________________________
Int_t Run()
{
  // start time (sec)
  gettimeofday(&fCurrentTimeStruct, NULL);
  fSecondsStart = fCurrentTimeStruct.tv_sec;

  // main loop
  while (!fgTerminationSignaled) {
    Int_t nSockets = 1;
    zmq_pollitem_t sockets[] = {
      { fZMQout, 0, ZMQ_POLLIN | ZMQ_POLLOUT, 0 },
    };

    Int_t rc = 0;
    errno = 0;

    Int_t outType = alizmq_socket_type(fZMQout);

    // poll sockets for data or conditions
    rc = zmq_poll(sockets, nSockets, fZMQtimeout); // poll sockets
    if (rc == -1 && errno == ETERM) {
      // this can only happen it the context was terminated, one of the sockets are
      // not valid or operation was interrupted
      Printf("zmq_poll was interrupted! rc = %i, %s", rc, zmq_strerror(errno));
      break;
    }

    if (alizmq_socket_type(fZMQout) == ZMQ_REP) {
      // request waiting
      if (sockets[0].revents & ZMQ_POLLIN) {
        printf("request in\n");
        aliZMQmsg message;
        alizmq_msg_recv(&message, fZMQout, 0);
        for (aliZMQmsg::iterator i = message.begin(); i != message.end(); ++i) {
          HandleRequest(i, fZMQout);
        }
        alizmq_msg_close(&message);
      }
    } // socket 0
    else if (sockets[0].revents & ZMQ_POLLOUT) {
      DoSend(fZMQout);
    }

  } // main loop

  // stop time (sec)
  gettimeofday(&fCurrentTimeStruct, NULL);
  fSecondsStop = fCurrentTimeStruct.tv_sec;

  {
    printf("number of messages sent    : %llu\n", fNumberOfMessagesSent);
    printf("out %lubytes, %.2f MB\n", fBytesOUT, (double)(fBytesOUT) / 1024. / 1024.);
  }

  return 0;
}

//_____________________________________________________________________
Int_t HandleControlSequence(aliZMQmsg::iterator block, void* socket)
{
  // handle control sequences in the request if any
  return 0;
}

//_____________________________________________________________________
Int_t HandleRequest(aliZMQmsg::iterator block, void* socket)
{
  printf("HandleRequest\n");
  HandleControlSequence(block, socket);
  return DoSend(socket);
}

//
struct NameHeader: public BaseDataTopic {
  static constexpr int nameSize{32};
  char _name[nameSize];
  NameHeader(char const* name, bool more = false):
    BaseDataTopic(sizeof(NameHeader), CharArr2uint64("NameHead"), CharArr2uint64("NONE"), more),
    _name{} {
      strncpy(_name, name, nameSize);
      //set padding length in last byte
      uint8_t* lastByte = reinterpret_cast<uint8_t*>(this) + sizeof(NameHeader) - 1;
      *lastByte = nameSize-strlen(name);
    };
};

struct Header {
  DataTopic dataHeader;
  NameHeader nameHeader;
  Header(char const* name, const char* id, const char* origin, uint64_t spec=0):
    dataHeader{id,origin,spec,CharArr2uint64("NONE"),true},
    nameHeader{name} {};
};

//______________________________________________________________________________
//this is a really stupid vector, push_back has no bounds checking!
//it tries to minimize reallocations and can release the underlying buffer
template<typename T>
class o2vec {
  public:
  o2vec(Header header, size_t capacity=0) :
    _buf(capacity?new T[capacity]:new T[1]), //don't want to deal with corner cases
    _capacity{capacity},
    _end(_buf.get()),
    _header(header)
    {}
  o2vec(const o2vec&) = delete;
  o2vec& operator=(const o2vec&) = delete;
  T* release() {
    return _buf.release();
  }
  size_t free() {
    return _buf.get()+_capacity-_end;
  }
  size_t size() {
    return _end-_buf.get();
  }
  size_t capacity() {
    return _capacity;
  }
  T& back() {
    return *(_end-1);
  }
  T* get() const {
    return _buf.get();
  }
  void push_back(T i) {
    *_end++=i;
  }
  void reserve(size_t elems) {
    if (elems>_capacity){
      T* tmp = new T[elems];
      memcpy(tmp,_buf.get(),size());
      _end=tmp+size();
      _buf.reset(tmp);
      _capacity=elems;
    };
  }
  Header* header() {return &_header;}
  void reserve_free(size_t elems) {
    if (elems>free()){
      reserve(size()+elems);
    };
  }
  private:
  std::unique_ptr<T[]> _buf;
  T* _end{nullptr};
  size_t _capacity{0};
  Header _header;
};

//______________________________________________________________________________
Int_t DoSend(void* socket)
{
  aliZMQmsg message;

  o2vec<int32_t> feventID(Header{"fEVID", "TRKPAREVID","ESD"});

  o2vec<float> fX(Header{"fX", "TRKPARX","ESD"});
  o2vec<float> fAlpha(Header{"fAlpha", "TRKPARAlpha","ESD"});
  o2vec<float> fY(Header{"fY", "TRKPARY","ESD"});
  o2vec<float> fZ(Header{"fZ", "TRKPARZ","ESD"});
  o2vec<float> fSnp(Header{"fSnp", "TRKPARSnp","ESD"});
  o2vec<float> fTgl(Header{"fTgl", "TRKPARTgl","ESD"});
  o2vec<float> fSigned1Pt(Header{"fSigned1Pt", "TRKPARSigned1Pt","ESD"});

  o2vec<float> fCYY(Header{"fCYY", "TRKCOVCYY","ESD"});
  o2vec<float> fCZY(Header{"fCZY", "TRKCOVCZY","ESD"});
  o2vec<float> fCZZ(Header{"fCZZ", "TRKCOVCZZ","ESD"});
  o2vec<float> fCSnpY(Header{"fCSnpY", "TRKCOVCSnpY","ESD"});
  o2vec<float> fCSnpZ(Header{"fCSnpZ", "TRKCOVCSnpZ","ESD"});
  o2vec<float> fCSnpSnp(Header{"fCSnpSnp", "TRKCOVCSnpSnp","ESD"});
  o2vec<float> fCTglY(Header{"fCTglY", "TRKCOVCTglY","ESD"});
  o2vec<float> fCTglZ(Header{"fCTglZ", "TRKCOVCTglZ","ESD"});
  o2vec<float> fCTglSnp(Header{"fCTglSnp", "TRKCOVCTglSnp","ESD"});
  o2vec<float> fCTglTgl(Header{"fCTglTgl", "TRKCOVCTglTgl","ESD"});
  o2vec<float> fC1PtY(Header{"fC1PtY", "TRKCOVC1PtY","ESD"});
  o2vec<float> fC1PtZ(Header{"fC1PtZ", "TRKCOVC1PtZ","ESD"});
  o2vec<float> fC1PtSnp(Header{"fC1PtSnp", "TRKCOVC1PtSnp","ESD"});
  o2vec<float> fC1PtTgl(Header{"fC1PtTgl", "TRKCOVC1PtTgl","ESD"});
  o2vec<float> fC1Pt21Pt2(Header{"fC1Pt21Pt2", "TRKCOVC1Pt21Pt2","ESD"});

  o2vec<float> fTPCinnerP(Header{"fTPCinnerP", "TRKEXTTPCinnerP","ESD"});
  o2vec<float> fFlags(Header{"fFlags", "TRKEXTFlags","ESD"});
  o2vec<float> fITSClusterMap(Header{"fITSClusterMap", "TRKEXTITSClsMap","ESD"});
  o2vec<float> fTPCncls(Header{"fTPCncls", "TRKEXTTPCncls","ESD"});
  o2vec<float> fTRDntracklets(Header{"fTRDntracklets", "TRKEXTTRDntrklts","ESD"});
  o2vec<float> fITSchi2Ncl(Header{"fITSchi2Ncl", "TRKEXTITSchi2Ncl","ESD"});
  o2vec<float> fTPCchi2Ncl(Header{"fTPCchi2Ncl", "TRKEXTTPCchi2Ncl","ESD"});
  o2vec<float> fTRDchi2(Header{"fTRDchi2", "TRKEXTTRDchi2","ESD"});
  o2vec<float> fTOFchi2(Header{"fTOFchi2", "TRKEXTTOFchi2","ESD"});
  o2vec<float> fTPCsignal(Header{"fTPCsignal", "TRKEXTTPCsignal","ESD"});
  o2vec<float> fTRDsignal(Header{"fTRDsignal", "TRKEXTTRDsignal","ESD"});
  o2vec<float> fTOFsignal(Header{"fTOFsignal", "TRKEXTTOFsignal","ESD"});
  o2vec<float> fLength(Header{"fLength", "TRKEXTLength","ESD"});

  for (std::string filename: fFilesIn) {

    if (fVerbose) {printf("opening file %si\n",filename.c_str());}
    TFile file(filename.c_str());
    if (file.IsOpen() == false) {
      printf("ERROR: cannot open %s",filename.c_str());
      continue;
    }

    TTree* esdTree = (TTree*)file.Get("esdTree");
    unique_ptr<AliESDEvent> esd(new AliESDEvent);
    esd->ReadFromTree(esdTree);
    size_t nev = esdTree->GetEntries();
    if (fVerbose) {printf("  %lu events\n",nev);}

    for (size_t iev = 0; iev < nev; ++iev) {
      esd->Reset();
      esdTree->GetEntry(iev);
      esd->ConnectTracks();

      // Tracks information
      size_t ntrk = esd->GetNumberOfTracks();

      if (fVerbose) {printf("  event: %lu tracks:%lu\n",iev,ntrk);}

      //naive preallocation scheme
      size_t elems = ntrk*nev*fFilesIn.size()+ntrk*nev;
      feventID.reserve(elems);
      fX.reserve(elems);
      fAlpha.reserve(elems);
      fY.reserve(elems);
      fZ.reserve(elems);
      fSnp.reserve(elems);
      fTgl.reserve(elems);
      fSigned1Pt.reserve(elems);
      fCYY.reserve(elems);
      fCZY.reserve(elems);
      fCZZ.reserve(elems);
      fCSnpY.reserve(elems);
      fCSnpZ.reserve(elems);
      fCSnpSnp.reserve(elems);
      fCTglY.reserve(elems);
      fCTglZ.reserve(elems);
      fCTglSnp.reserve(elems);
      fCTglTgl.reserve(elems);
      fC1PtY.reserve(elems);
      fC1PtZ.reserve(elems);
      fC1PtSnp.reserve(elems);
      fC1PtTgl.reserve(elems);
      fC1Pt21Pt2.reserve(elems);
      fTPCinnerP.reserve(elems);
      fFlags.reserve(elems);
      fITSClusterMap.reserve(elems);
      fTPCncls.reserve(elems);
      fTRDntracklets.reserve(elems);
      fITSchi2Ncl.reserve(elems);
      fTPCchi2Ncl.reserve(elems);
      fTRDchi2.reserve(elems);
      fTOFchi2.reserve(elems);
      fTPCsignal.reserve(elems);
      fTRDsignal.reserve(elems);
      fTOFsignal.reserve(elems);
      fLength.reserve(elems);

      for (size_t itrk = 0; itrk < ntrk; ++itrk) {
        AliESDtrack* track = esd->GetTrack(itrk);
        track->SetESDEvent(esd.get());

        feventID.push_back(iev);
        fX.push_back(track->GetX());
        fAlpha.push_back(track->GetAlpha());
        fY.push_back(track->GetY());
        fZ.push_back(track->GetZ());
        fSnp.push_back(track->GetSnp());
        fTgl.push_back(track->GetTgl());
        fSigned1Pt.push_back(track->GetSigned1Pt());

        fCYY.push_back(track->GetSigmaY2());
        fCZY.push_back(track->GetSigmaZY());
        fCZZ.push_back(track->GetSigmaZ2());
        fCSnpY.push_back(track->GetSigmaSnpY());
        fCSnpZ.push_back(track->GetSigmaSnpZ());
        fCSnpSnp.push_back(track->GetSigmaSnp2());
        fCTglY.push_back(track->GetSigmaTglY());
        fCTglZ.push_back(track->GetSigmaTglZ());
        fCTglSnp.push_back(track->GetSigmaTglSnp());
        fCTglTgl.push_back(track->GetSigmaTgl2());
        fC1PtY.push_back(track->GetSigma1PtY());
        fC1PtZ.push_back(track->GetSigma1PtZ());
        fC1PtSnp.push_back(track->GetSigma1PtSnp());
        fC1PtTgl.push_back(track->GetSigma1PtTgl());
        fC1Pt21Pt2.push_back(track->GetSigma1Pt2());

        const AliExternalTrackParam* intp = track->GetTPCInnerParam();

        fTPCinnerP.push_back(intp ? intp->GetP() : 0);
        fFlags.push_back(track->GetStatus());
        fITSClusterMap.push_back(track->GetITSClusterMap());
        fTPCncls.push_back(track->GetTPCNcls());
        fTRDntracklets.push_back(track->GetTRDntracklets());
        fITSchi2Ncl.push_back(track->GetITSNcls() ? track->GetITSchi2() / track->GetITSNcls() : 0);
        fTPCchi2Ncl.push_back(track->GetTPCNcls() ? track->GetTPCchi2() / track->GetTPCNcls() : 0);
        fTRDchi2.push_back(track->GetTRDchi2());
        fTOFchi2.push_back(track->GetTOFchi2());
        fTPCsignal.push_back(track->GetTPCsignal());
        fTRDsignal.push_back(track->GetTRDsignal());
        fTOFsignal.push_back(track->GetTOFsignal());
        fLength.push_back(track->GetIntegratedLength());

      }//track loop
    }//event loop
  }//file loop

  int rc{0};
  rc = alizmq_msg_add_array(&message, std::move(feventID));
  rc = alizmq_msg_add_array(&message, std::move(fX));
  rc = alizmq_msg_add_array(&message, std::move(fAlpha));
  rc = alizmq_msg_add_array(&message, std::move(fY));
  rc = alizmq_msg_add_array(&message, std::move(fZ));
  rc = alizmq_msg_add_array(&message, std::move(fSnp));
  rc = alizmq_msg_add_array(&message, std::move(fTgl));
  rc = alizmq_msg_add_array(&message, std::move(fSigned1Pt));
  rc = alizmq_msg_add_array(&message, std::move(fCYY));
  rc = alizmq_msg_add_array(&message, std::move(fCZY));
  rc = alizmq_msg_add_array(&message, std::move(fCZZ));
  rc = alizmq_msg_add_array(&message, std::move(fCSnpY));
  rc = alizmq_msg_add_array(&message, std::move(fCSnpZ));
  rc = alizmq_msg_add_array(&message, std::move(fCSnpSnp));
  rc = alizmq_msg_add_array(&message, std::move(fCTglY));
  rc = alizmq_msg_add_array(&message, std::move(fCTglZ));
  rc = alizmq_msg_add_array(&message, std::move(fCTglSnp));
  rc = alizmq_msg_add_array(&message, std::move(fCTglTgl));
  rc = alizmq_msg_add_array(&message, std::move(fC1PtY));
  rc = alizmq_msg_add_array(&message, std::move(fC1PtZ));
  rc = alizmq_msg_add_array(&message, std::move(fC1PtSnp));
  rc = alizmq_msg_add_array(&message, std::move(fC1PtTgl));
  rc = alizmq_msg_add_array(&message, std::move(fC1Pt21Pt2));
  rc = alizmq_msg_add_array(&message, std::move(fTPCinnerP));
  rc = alizmq_msg_add_array(&message, std::move(fFlags));
  rc = alizmq_msg_add_array(&message, std::move(fITSClusterMap));
  rc = alizmq_msg_add_array(&message, std::move(fTPCncls));
  rc = alizmq_msg_add_array(&message, std::move(fTRDntracklets));
  rc = alizmq_msg_add_array(&message, std::move(fITSchi2Ncl));
  rc = alizmq_msg_add_array(&message, std::move(fTPCchi2Ncl));
  rc = alizmq_msg_add_array(&message, std::move(fTRDchi2));
  rc = alizmq_msg_add_array(&message, std::move(fTOFchi2));
  rc = alizmq_msg_add_array(&message, std::move(fTPCsignal));
  rc = alizmq_msg_add_array(&message, std::move(fTRDsignal));
  rc = alizmq_msg_add_array(&message, std::move(fTOFsignal));
  rc = alizmq_msg_add_array(&message, std::move(fLength));

  fBytesOUT = alizmq_msg_send(&message, socket, 0);
  if (fBytesOUT<=0) {
    printf("ERROR sending: %s\n",zmq_strerror(zmq_errno()));
  }
  fNumberOfMessagesSent++;
  SetLastPushBackTime();
  alizmq_msg_close(&message);

  if (fOnce) { fgTerminationSignaled=true; }

  return fBytesOUT;
}

//______________________________________________________________________________
void SetLastPushBackTime()
{
  gettimeofday(&fCurrentTimeStruct, NULL);
  fLastPushBackTime = 1000 * fCurrentTimeStruct.tv_sec + fCurrentTimeStruct.tv_usec / 1000;
}

//_______________________________________________________________________________________
Int_t InitZMQ()
{
  // init or reinit stuff
  Int_t rc = 0;
  rc = alizmq_socket_init(fZMQout, fZMQcontext, fZMQconfigOUT.Data(), fZMQtimeout, fZMQmaxQueueSize, fZMQlingerTime);
  printf("out: (%s) %s\n", alizmq_socket_name(rc), fZMQconfigOUT.Data());
  return 0;
}

//______________________________________________________________________________
Int_t ProcessOptionString(TString arguments, Bool_t verbose)
{
  // process passed options
  Int_t nOptions = 0;
  std::unique_ptr<aliStringVec> options{ AliOptionParser::TokenizeOptionString(arguments) };
  for (auto i : *options) {
    const TString& option = i.first;
    const TString& value = i.second;
    if (option.EqualTo("ZMQconfigOUT") || option.EqualTo("out")) {
      fZMQconfigOUT = value;
      fReconfigureZMQ = true;
    }
    else if (option.EqualTo("in")) {
        Ssiz_t from{0};
        TString token{};
        while (value.Tokenize(token,from,",")) {
          fFilesIn.emplace_back(token.Data());
        }
    }
    else if (option.EqualTo("Verbose")) {
      fVerbose = value.EqualTo("0") ? kFALSE : kTRUE;
    }
    else if (option.EqualTo("ZMQmaxQueueSize")) {
      fZMQmaxQueueSize = value.Atoi();
      fReconfigureZMQ = true;
    }
    else if (option.EqualTo("ZMQtimeout")) {
      fZMQtimeout = value.Atoi();
      fReconfigureZMQ = true;
    }
    else if (option.EqualTo("id")) {
      fID = value.Data();
    }
    else if (option.EqualTo("once")) {
      fOnce = value.Atoi();
    }
    else if (option.EqualTo("ZMQlinger")) {
      fZMQlingerTime = value.Atoi();
    }
    else {
      Printf("unrecognized option |%s|", option.Data());
      nOptions = -1;
      break;
    }
    nOptions++;
  }

  if (nOptions < 1)
    fReconfigureZMQ = false;
  if (fReconfigureZMQ && (InitZMQ() < 0)) {
    Printf("failed ZMQ init");
    return -1;
  }
  fReconfigureZMQ = false;

  if (fRequestTimeout < 100)
    printf("WARNING: setting the socket timeout to %lu ms can be dagerous,\n"
           "         choose something more realistic or leave the default as it is\n",
           fRequestTimeout);

  return nOptions;
}

//_______________________________________________________________________________________
int main(Int_t argc, char** argv)
{
  Int_t mainReturnCode = 0;

  // catch signals
  if (signal(SIGHUP, sig_handler) == SIG_ERR)
    printf("\ncan't catch SIGHUP\n");
  if (signal(SIGINT, sig_handler) == SIG_ERR)
    printf("\ncan't catch SIGINT\n");
  if (signal(SIGQUIT, sig_handler) == SIG_ERR)
    printf("\ncan't catch SIGQUIT\n");
  if (signal(SIGTERM, sig_handler) == SIG_ERR)
    printf("\ncan't catch SIGTERM\n");

  // the context
  fZMQcontext = alizmq_context();
  if (!fZMQcontext) {
    printf("could not init the ZMQ context\n");
    return 1;
  }

  // process args
  TString argString = AliOptionParser::GetFullArgString(argc, argv);
  if (ProcessOptionString(argString, kTRUE) <= 0) {
    printf("%s", fUSAGE);
    return 1;
  }

  Run();

  // destroy ZMQ sockets
  zmq_close(fZMQout);
  zmq_ctx_destroy(fZMQcontext);
  return mainReturnCode;
}

//_______________________________________________________________________________________
template <typename T>
void alizmq_deleteArray(void* data, void*) {
  delete [] static_cast<T*>(data);
}

template <typename T>
int alizmq_msg_add_array(aliZMQmsg* message, o2vec<T>&& vec)
{
  //add a frame to the mesage
  int rc = 0;

  size_t nelems = vec.size();
  T* array = vec.release();
  Header* topic = vec.header();

  //prepare topic msg
  zmq_msg_t* topicMsg = new zmq_msg_t;
  rc = zmq_msg_init_size( topicMsg, sizeof(*topic));
  if (rc<0) {
    zmq_msg_close(topicMsg);
    delete topicMsg;
    return -1;
  }
  memcpy(zmq_msg_data(topicMsg),topic,sizeof(*topic));

  size_t size = nelems*sizeof(T);

  //prepare data msg
  zmq_msg_t* dataMsg = new zmq_msg_t;
  rc = zmq_msg_init_data( dataMsg, array, size, alizmq_deleteArray<T>, nullptr);
  if (rc<0) {
    zmq_msg_close(topicMsg);
    zmq_msg_close(dataMsg);
    delete topicMsg;
    delete dataMsg;
    return -1;
  }

  static_cast<DataTopic*>(zmq_msg_data(topicMsg))->SetPayloadSize(size);

  //add the frame to the message
  message->push_back(std::make_pair(topicMsg,dataMsg));
  return message->size();
}

