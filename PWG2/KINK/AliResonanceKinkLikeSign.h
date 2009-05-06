#ifndef ALIRESONANCEKINKLIKESIGN_H
#define ALIRESONANCEKINKLIKESIGN_H

/*  See cxx source for full Copyright notice */

//--------------------------------------------------------------------------------
//                   class AliResonanceKinkLikeSign
//         This task is an example of an analysis task
//        for producing a like-sign background for resonances having at least one 
//        kaon-kink in their decay products. 
//        Author: Paraskevi Ganoti, University of Athens (pganoti@phys.uoa.gr)
//---------------------------------------------------------------------------------
class TF1;
class TTree;
class AliESDEvent;
class AliPID;
class AliAnalysisTaskSE;

class AliResonanceKinkLikeSign : public AliAnalysisTaskSE {
 public:
  AliResonanceKinkLikeSign();
  AliResonanceKinkLikeSign(const char *name);
  virtual ~AliResonanceKinkLikeSign() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  Float_t GetSigmaToVertex(AliESDtrack* esdTrack) const ; 
  const AliESDVertex *GetEventVertex(const AliESDEvent* esd) const;
  
 private:
  AliESDEvent *fESD;    //! ESD object
  TList       *fListOfHistos; //! List 
  TF1         *f1;
  TF1         *f2;  
  TH1D        *fNegKaonLikeSign; //! negative spectrum
  
  AliResonanceKinkLikeSign(const AliResonanceKinkLikeSign&); // not implemented
  AliResonanceKinkLikeSign& operator=(const AliResonanceKinkLikeSign&); // not implemented

  ClassDef(AliResonanceKinkLikeSign, 1); // example of analysis
};

#endif
