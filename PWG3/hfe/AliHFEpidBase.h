#ifndef ALIHFEPIDBASE_H
#define ALIHFEPIDBASE_H
 
/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

/* $Id$ */ 

//
// Base Class for Detector PID Objects
// For more information see the implementation file
//

 #ifndef ROOT_TNamed
 #include <TNamed.h>
 #endif

class TList;
class AliAODpidUtil;
class AliESDpid;
class AliVParticle;
class AliMCParticle;
class AliHFEpidQAmanager;

class AliHFEpidObject{
  public:
    typedef enum{ 
      kESDanalysis,
      kAODanalysis
    }AnalysisType_t;
    AliHFEpidObject():
      fkRecTrack(NULL), 
      fAnalysisType(kESDanalysis),
      fAbInitioPID(-1),
      fCentrality(99.)
      {
      }
    AliHFEpidObject(const AliHFEpidObject &ref):
      fkRecTrack(ref.fkRecTrack), 
      fAnalysisType(ref.fAnalysisType),
      fAbInitioPID(ref.fAbInitioPID),
      fCentrality(ref.fCentrality)
      {
      }
    AliHFEpidObject &operator=(const AliHFEpidObject &ref);
    ~AliHFEpidObject(){};

    void SetRecTrack(const AliVParticle * recTrack) {fkRecTrack = recTrack; }
    void SetMCTrack(const AliVParticle * mcTrack);
    void SetAnalysisType(AnalysisType_t type) { fAnalysisType = type; }
    void SetAbInitioPID(Int_t abInitioPID) { fAbInitioPID = abInitioPID; }
    void SetCentrality(Float_t centrality) { fCentrality = centrality; }

    const AliVParticle *GetRecTrack() const { return fkRecTrack; }
    Int_t GetAbInitioPID() const { return fAbInitioPID; }
    Float_t GetCentrality() const { return fCentrality; }
    Bool_t IsAODanalysis() const { return fAnalysisType == static_cast<UChar_t>(kAODanalysis); }
    Bool_t IsESDanalysis() const { return fAnalysisType == static_cast<UChar_t>(kESDanalysis); }

  private:
    const AliVParticle *fkRecTrack;    // Reconstructed track
    UChar_t fAnalysisType;      // Analysis Mode (ESD or AOD)
    Int_t fAbInitioPID;         // AbInitio PID
    Float_t fCentrality;        // Centrality Information
};

class AliHFEpidBase : public TNamed{
  public:
    AliHFEpidBase();
    AliHFEpidBase(const Char_t *name);
    AliHFEpidBase(const AliHFEpidBase &c);
    AliHFEpidBase &operator=(const AliHFEpidBase &c);
    virtual ~AliHFEpidBase() {};
    // Framework functions that have to be implemented by the detector PID classes
    virtual Bool_t InitializePID() = 0;
    virtual Int_t IsSelected(const AliHFEpidObject *track, AliHFEpidQAmanager *pidqa = NULL) const = 0;

    Bool_t HasMCData() const { return TestBit(kHasMCData); };

    void SetESDpid(AliESDpid * const pid) { fESDpid = pid; }
    void SetAODpid(AliAODpidUtil * const pid) { fAODpid = pid; }
    void SetHasMCData(Bool_t hasMCdata = kTRUE) { SetBit(kHasMCData,hasMCdata); };

    AliESDpid *GetESDpid() const { return fESDpid; }; 

  protected:
    AliESDpid *fESDpid;                         //! ESD PID object
    AliAODpidUtil *fAODpid;                     //! AOD PID object
    void Copy(TObject &ref) const;

  private:
    enum{
      kHasMCData = BIT(14)
    };

    ClassDef(AliHFEpidBase, 2)      // Base class for detector Electron ID
};
#endif
