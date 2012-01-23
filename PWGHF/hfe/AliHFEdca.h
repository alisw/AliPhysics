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
//
// Class for checking impact parameter (DCA) study 
// + study DCA in rphi (xy) and z
// + resolution and pull
// + handle both MC and data 
// + add plugin for primary vertex 
//

#ifndef ALIHFEDCA_H
#define ALIHFEDCA_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

class TChain;
class TTree;
class TFile;

class TString;
class TList;

class TObjArray;
class AliStack;
class AliMCEvent;
class AliMCVertex;

class AliESDEvent;
class AliESDtrack;
class AliESDVertex;

class AliHFEdca : public TObject{

 public:  
  enum{
    kPDGelectron = 11,
    kPDGmuon = 13,
    kPDGpion = 211,
    kPDGkaon = 321,
    kPDGproton = 2212
  };
 
  enum{
    kNParticles = 12,
    kNPtBins = 50,   
    kNDcaVar = 2, 
    kNVertexVar = 3,  
    kNPullVar = 2
  };

  AliHFEdca(); // default constructor
  AliHFEdca(const AliHFEdca &ref); // copy constructor
  AliHFEdca &operator=(const AliHFEdca &ref); // assignment operator
  virtual ~AliHFEdca(); // destructor

  void Initialize();
  void CreateHistogramsPull(TList *pullList);  
  void CreateHistogramsResidual(TList *residualList);  
  void CreateHistogramsDca(TList *dcaList);  

  void CreateHistogramsKfDca(TList *kfDcaList);  

  void CreateHistogramsDataDca(TList *dataDcaList);  
  void CreateHistogramsDataPull(TList *dataPullList);  

  void CreateHistogramsVertex(TList *vertexList);  
  void CreateHistogramsDataVertex(TList *vertexList);  

  void CreateHistogramsPid(TList *pidList);
  void CreateHistogramsDataPid(TList *pidList);

  void CreateHistogramsHfeDca(TList *hfeDcaList);
  void CreateHistogramsHfeDataDca(TList *hfeDataDcaList);

  
  void InitAnalysis()const;  
  void FillHistogramsDca(const AliESDEvent * const esdEvent,  const AliESDtrack *const track,  AliMCEvent *const mcEvent);
  void FillHistogramsVtx(const AliESDEvent * const esdEvent,  const AliMCEvent *const mcEvent);
  void FillHistogramsPid(const AliESDtrack *track, const  AliMCEvent * const mcEvent);

  void FillHistogramsKfDca(const AliESDEvent * const esdEvent,  const AliESDtrack *const track,  const AliMCEvent *const mcEvent);

  void FillHistogramsDataDca(const AliESDEvent * const esdEvent,  const AliESDtrack * const track, const AliESDVertex * const vtxESDSkip);
  void FillHistogramsDataVtx(const AliESDEvent * const esdEvent);
  void FillHistogramsDataPid(const AliESDtrack * const track);

  void FillHistogramsHfeDca(const AliESDEvent *const esdEvent,  const AliESDtrack * const track,  const AliMCEvent * const mcEvent);
  void FillHistogramsHfeDataDca(const AliESDEvent * const esdEvent, const AliESDtrack * const track, const AliESDVertex * const vtxESDSkip);


  void ApplyExtraCuts(const AliESDEvent * const esdEvent, Int_t nMinPrimVtxContributor);

  void PostAnalysis() const;

  Int_t GetCombinedPid(const AliESDtrack * const track);

 private:   


  static const Char_t* fgkParticles[kNParticles];  // particle names
  static const Int_t fgkPdgParticle[kNParticles-2]; // identified particle's name
  static const Int_t fgkColorPart[kNParticles]; // colors for particles

  static const Float_t fgkPtIntv[kNPtBins+1];  // pt intervals

  static const Char_t* fgkDcaVar[kNDcaVar];  // dca variables
  static const Char_t* fgkDcaVarTitle[kNDcaVar]; // titles for dca variables

  static const Char_t* fgkVertexVar[kNVertexVar];  // dca variables
  static const Char_t* fgkVertexVarTitle[kNVertexVar]; // titles for dca variables

  static const Char_t* fgkResDcaVar[kNDcaVar];  // dca variables
  static const Char_t* fgkResDcaVarTitle[kNDcaVar]; // titles for dca variables

  static const Char_t* fgkPullDcaVar[kNPullVar];  // pull variables
  static const Char_t* fgkPullDcaVarTitle[kNPullVar]; // titles for pull variables
  static const Char_t* fgkPullDataDcaVarTitle[kNPullVar]; // titles for pull variables

  TH1F* fHistDcaXYRes[kNParticles][kNPtBins];  //! residuals in XY
  TH1F* fHistDcaZRes[kNParticles][kNPtBins];   //! residuals in Z

  TH1F* fHistDcaXYPull[kNParticles][kNPtBins]; //! pulls XY
  TH1F* fHistDcaZPull[kNParticles][kNPtBins];  //! pulls Z

  TH1F* fHistDcaXY[kNParticles][kNPtBins]; //! dca XY
  TH1F* fHistDcaZ[kNParticles][kNPtBins];  //! dca Z

  TH1F* fHistEPDcaXYRes[kNParticles-2][kNPtBins];  //! residuals in XY with esd pid
  TH1F* fHistEPDcaZRes[kNParticles-2][kNPtBins];   //! residuals in Z with esd pid

  TH1F* fHistEPDcaXYPull[kNParticles-2][kNPtBins]; //! pulls XY with esd pid
  TH1F* fHistEPDcaZPull[kNParticles-2][kNPtBins];  //! pulls Z with esd pid

  TH1F* fHistEPDcaXY[kNParticles-2][kNPtBins]; //! dca XY with esd pid
  TH1F* fHistEPDcaZ[kNParticles-2][kNPtBins];  //! dca Z with esd pid

  TH1F* fHistKFDcaXY[kNParticles][kNPtBins]; //! KF dca XY
  TH1F* fHistKFDcaZ[kNParticles][kNPtBins];  //! KF dca Z
  
  TH1F* fHistDataDcaXY[kNParticles][kNPtBins]; //! data dca XY
  TH1F* fHistDataDcaZ[kNParticles][kNPtBins];  //! data dca Z
  TH1F* fHistDataWoDcaXY[kNParticles][kNPtBins]; //! data dca XY w/o current trk
  TH1F* fHistDataWoDcaZ[kNParticles][kNPtBins];  //! data dca Z w/o current trk

  TH1F* fHistDataDcaXYPull[kNParticles][kNPtBins]; //! data pull dca XY
  TH1F* fHistDataDcaZPull[kNParticles][kNPtBins];  //! data pull dca Z 
  TH1F* fHistDataWoDcaXYPull[kNParticles][kNPtBins]; //! data pull dca XY w/o current trk
  TH1F* fHistDataWoDcaZPull[kNParticles][kNPtBins];  //! data pull dca Z  w/o current trk

  TH1F* fHistMCvertex[kNVertexVar];    //! vertex MC
  TH1F* fHistESDvertex[kNVertexVar];   //! vertex ESD
  TH1F* fHistDatavertex[kNVertexVar];  //! vertex Data
  
  TH1F* fHistMcPid[kNParticles];      //! MC pid pt spectra
  TH1F* fHistEsdPid[kNParticles];     //! ESD pid pt spectra

  TH1F *fHistDataEsdPid[kNParticles];    //! Data ESD pid

  // HFE pid part
  // MC
  TH1F* fHistHPDcaXYRes[2][kNPtBins];  //! residuals in XY
  TH1F* fHistHPDcaZRes[2][kNPtBins];   //! residuals in Z
  TH1F* fHistHPDcaXYPull[2][kNPtBins]; //! pulls XY
  TH1F* fHistHPDcaZPull[2][kNPtBins];  //! pulls Z
  TH1F* fHistHPDcaXY[2][kNPtBins];     //! dca XY
  TH1F* fHistHPDcaZ[2][kNPtBins];      //! dca Z

  TH1F* fHistHfePid[2][2];             // ! HFE pid pt spectra only for electrons

  // Data
  TH1F* fHistHPDataDcaXY[2][kNPtBins];     //! data dca XY with HFE pid
  TH1F* fHistHPDataDcaZ[2][kNPtBins];      //! data dca Z with HFE pid
  TH1F* fHistHPDataDcaXYPull[2][kNPtBins]; //! data pull dca XY
  TH1F* fHistHPDataDcaZPull[2][kNPtBins];  //! data pull dca Z 

  TH1F *fHistDataHfePid[2];            //! Data HFE pid

  TH1I* fStat;                         //! counting diff of dca calculated from HF particle and ESD
  ClassDef(AliHFEdca, 1);
};

#endif
