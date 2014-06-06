#ifndef ALIVVEVENTFRIEND_H
#define ALIVVEVENTFRIEND_H

class AliVVVZEROfriend;
class AliVVTZEROfriend;

//_____________________________________________________________________________
class AliVVfriend {
public:
  AliVVfriend();
  virtual ~AliVVfriendEvent();

  virtual Int_t GetNumberOfTracks() const {return 0;}
  AliVVfriendTrack* GetTrack(Int_t /*i*/) const {return NULL;}
  virtual Int_t GetEntriesInTracks() const {return 0;}
  
  virtual AliESDVZEROfriend *GetVZEROfriend(){ return NULL; }
  virtual AliESDTZEROfriend *GetTZEROfriend(){ return NULL; }

  virtual void Ls() {}

  virtual void Reset() {}
  
  Bool_t TestSkipBit() {return kFALSE;}

  virtual Int_t GetNclustersTPC(UInt_t /*sector*/) const {return 0;}
  virtual Int_t GetNclustersTPCused(UInt_t /*sector*/) const {return 0;}
  
  //virtual void AddTrack(const AliVVfriendTrack *t) {}
  //virtual void AddTrackAt(const AliVVfriendTrack* /*t*/, Int_t /*i*/) {}
  //virtual void SetVZEROfriend(AliESDVZEROfriend* /*obj*/) {}
  //virtual void SetTZEROfriend(AliESDTZEROfriend * obj) {}
  //void SetSkipBit(Bool_t skip){}

private: 
  AliVVfriend(const AliVVfriend &);
  AliVVfriend& operator=(const AliVVfriend& esd);  
};

#endif


