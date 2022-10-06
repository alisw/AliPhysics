/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */
#ifndef AliMCEventContainer_H
#define AliMCEventContainer_H

#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"


class AliMCEventContainer : public TObject
{
public:
    AliMCEventContainer(); // default constructor
    AliMCEventContainer(const AliMCEventContainer &source); // copy constructor
    virtual ~AliMCEventContainer();//destructor.
    
    void Fill(AliMultSelection* MultSelection, AliAnalysisUtils* Utils, AliAODEvent* event);
    void Reset();
    Bool_t SelectVertex2015pp(AliAODEvent *aod,
                                Bool_t checkSPDres = kTRUE, //enable check on vtx resolution 
                                Bool_t *SPDandTrkExists = NULL, //ask for both trk and SPD vertex
                                Bool_t *checkProximity = NULL); //apply cut on relative position of spd and trk verteces
    
    Bool_t IsGoodSPDvertexRes(const AliAODVertex * spdVertex = NULL);
	
	Bool_t GetVertexStatus(const AliAODVertex *vertex);
    
    //setter 
    void SetPrimaryVertexPos(Double_t Pos[3])
    {
        fPrimVertex[0] = Pos[0];
        fPrimVertex[1] = Pos[1];
        fPrimVertex[2] = Pos[2];
    }
    
    void SetPrimaryVertexCov(Double_t Cov[6])
    {
        fPrimVertexCov[0] = Cov[0];
        fPrimVertexCov[1] = Cov[1];
        fPrimVertexCov[2] = Cov[2];
        fPrimVertexCov[3] = Cov[3];
        fPrimVertexCov[4] = Cov[4];
        fPrimVertexCov[5] = Cov[5];
    }
    
    void SetPrimaryVertexMCPos(Double_t Pos[3])
    {
        fPrimVertexMC[0] = Pos[0];
        fPrimVertexMC[1] = Pos[1];
        fPrimVertexMC[2] = Pos[2];
    }
    
    void SetTriggerMask(UInt_t mask)    { fTriggerMask = mask ;}
    
    void SetIsCollisionCandidate    (Bool_t value)      { fIsCollisionCandidate = value ; }
    void SetEventNumber             (ULong64_t value)   { fEventNumber          = value ; }
    
    void SetMCMultiplicity          (Int_t  value)      { fMCMultiplicity       = value ; }
    
    
    //getter
    Bool_t      IsCollisionCandidate        ()           const { return fIsCollisionCandidate;       }
    UInt_t      GetTriggerMask              ()           const { return fTriggerMask;                }
    
    Double_t    GetMagneticField            ()           const { return fMagField;                   }
    Int_t       GetRunNumber                ()           const { return fRunNumber;                  }
    ULong64_t   GetEventNumber              ()           const { return fEventNumber;                }
    Double_t    GetPrimVertex               (Int_t i)    const { return fPrimVertex[i];              }
    Double_t    GetPrimVertexCov            (Int_t i)    const { return fPrimVertexCov[i];           }
    Double_t    GetPrimVertexMC             (Int_t i)    const { return fPrimVertexMC[i];            }
    
    Bool_t      HasVertex                   ()           const { return fHasVertex;                  }
    Bool_t      GetProximityCut             ()           const { return fProximityCut;               }
    Bool_t      HasSPDANDTrkVtx             ()           const { return fHasSPDANDTrkVtx;            }
    Bool_t      IsVertexSelected2015pp      ()           const { return fVertexSelected2015pp;       }
    
    Float_t     GetV0MPercCorr              ()           const { return fV0MPercCorr;                }
    Float_t     GetV0MPercentile            ()           const { return fV0MPercentile;              }
    Float_t     GetV0APercentile            ()           const { return fV0APercentile;              }
    Float_t     GetV0CPercentile            ()           const { return fV0CPercentile;              }
    
    Float_t     GetV0ASignal                ()           const { return fV0ASignal;                  }
    Float_t     GetV0CSignal                ()           const { return fV0CSignal;                  }
    
    Float_t     GetSPDTrackletsPerc         ()           const { return fSPDTrackletsPerc;           }
    Float_t     GetSPDTrackletsInEta08Perc  ()           const { return fSPDTrackletsInEta08Perc;    }
    Float_t     GetSPDTrackletsInEta05Perc  ()           const { return fSPDTrackletsInEta05Perc;    }
    
    Float_t     GetSPDTracklets             ()           const { return fSPDTracklets;               }
    Float_t     GetSPDTrackletsInEta08      ()           const { return fSPDTrackletsInEta08;        }
    Float_t     GetSPDTrackletsInEta05      ()           const { return fSPDTrackletsInEta05;        }
    
    Int_t       GetMCMultiplicity           ()           const { return fMCMultiplicity;             }
    
    Bool_t      IsSPDClusterVsTrackletBG    ()           const { return fIsSPDClusterVsTrackletBG;   }
    Bool_t      IsIncompleteDAQ             ()           const { return fIsIncompleteDAQ;            }
    Bool_t      IsINELgtZERO                ()           const { return fIsINELgtZERO;               }
    Bool_t      HasNoInconsistentVertices   ()           const { return fHasNoInconsistentVertices;  }
    Bool_t      IsNotAsymmetricInVZERO      ()           const { return fIsNotAsymmetricInVZERO;     }
    
    Bool_t      IsPileupMV                  ()           const { return fIsPileupMV;                 }
    Bool_t      IsPileupOOB                 ()           const { return fIsPileupOOB;                }
    Bool_t      IsPileupFromSPD             ()           const { return fIsPileupFromSPD;            }
    Bool_t      IsPileupFromSPDInMultBins   ()           const { return fIsPileupFromSPDInMultBins;  }
    Bool_t      IsPileupInMultBins          ()           const { return fIsPileupInMultBins;         }
    
private:
    //----------------GENERAL-INFO--------------------
    Bool_t      fIsCollisionCandidate;      // is collision candidate
    UInt_t      fTriggerMask;               // save full info for checking later
    
    Double_t    fMagField;                  //
    Int_t       fRunNumber;                 //
    ULong64_t   fEventNumber;               //
    
    Char_t      fPrimVertexType;            //
    Double_t    fPrimVertex     [3];        // [0] = X ; [1] = Y ; [2] = Z
    Double_t    fPrimVertexCov  [6];        // [0] = XX, [1] = XY, [2] = XZ ; [3] = YY, [4] = YZ, [5] = ZZ
    Double_t    fPrimVertexMC   [3];        // [0] = X ; [1] = Y ; [2] = Z
    Bool_t      fHasVertex;                 //
    Bool_t      fProximityCut;              //
    Bool_t      fHasSPDANDTrkVtx;           //
    Bool_t      fVertexSelected2015pp;      //
    //--------------MULTIPLICITY-INFO-----------------
    //------------------------------------------------
    // VZERO info
    Float_t     fV0MPercCorr;               // VZERO Corrected Percentile   - AliMultSelection
    Float_t     fV0MPercentile;             // VZERO Percentile             - AliMultSelection
    Float_t     fV0APercentile;             // VZERO Percentile             - AliMultSelection
    Float_t     fV0CPercentile;             // VZERO Percentile             - AliMultSelection
    
    Float_t     fV0ASignal;                 //
    Float_t     fV0CSignal;                 //
    //------------------------------------------------
    //SPD info
    Float_t     fSPDTrackletsPerc;          // SPD Multiplicity             - SPD Tracklets
    Float_t     fSPDTrackletsInEta08Perc;   // SPD Multiplicity             - SPD Tracklets
    Float_t     fSPDTrackletsInEta05Perc;   // SPD Multiplicity             - SPD Tracklets
    Float_t     fSPDTracklets;              // SPD Multiplicity             - SPD Tracklets
    Float_t     fSPDTrackletsInEta08;       // SPD Multiplicity             - SPD Tracklets
    Float_t     fSPDTrackletsInEta05;       // SPD Multiplicity             - SPD Tracklets
    //------------------------------------------------
    //MC info
    Int_t       fMCMultiplicity;            // 
    //---------------PILE-UP-INFO---------------------
    Bool_t      fIsSPDClusterVsTrackletBG;  //
    Bool_t      fIsIncompleteDAQ;           //
    Bool_t      fIsINELgtZERO;              //
    Bool_t      fHasNoInconsistentVertices; //
    Bool_t      fIsNotAsymmetricInVZERO;    //
    Bool_t		fHasGoodVertex2016;			//
    //---------------PILE-UP-INFO---------------------
    Bool_t      fIsPileupMV;                //
    Bool_t      fIsPileupOOB;               //
    Bool_t      fIsPileupFromSPD;           //
    Bool_t      fIsPileupFromSPDInMultBins; //
    Bool_t      fIsPileupInMultBins;        //
    //------------------------------------------------
    ClassDef(AliMCEventContainer, 1);
};
#endif
