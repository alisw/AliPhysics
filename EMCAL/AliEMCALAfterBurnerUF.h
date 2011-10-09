#ifndef ALIEMCALAFTERBURNERUF_H
#define ALIEMCALAFTERBURNERUF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
//  After-burner for the EMCAL cluster unfolding algorithm
//
//  See cxx for details on how to use it
//
//  Author: Olga Driga (SUBATECH)
//

// --- ROOT system ---
class TObjArray;
class TClonesArray;

// --- Standard library ---

// --- AliRoot header files ---
class AliEMCALGeometry;
class AliEMCALUnfolding;
class AliVCaloCells;

class AliEMCALAfterBurnerUF {

  public:
    AliEMCALAfterBurnerUF();
    AliEMCALAfterBurnerUF(Float_t logWeight, Float_t locMaxCut, Float_t minEcut);
    virtual ~AliEMCALAfterBurnerUF();

  private:
    AliEMCALAfterBurnerUF(const AliEMCALAfterBurnerUF & uf) ; // cpy ctor not needed, put here to avoid compilation warning 
    AliEMCALAfterBurnerUF & operator = (const AliEMCALAfterBurnerUF & uf) ;//cpy assignment, put here to avoid compilation warning 
  
  public:
    virtual void Clear();
    virtual void Init();
    virtual void RecPoints2Clusters(TObjArray *clusArray);
    virtual void UnfoldClusters(TObjArray *clusArray, AliVCaloCells *cellsEMCAL);  // does the job
  
    // getters and setters
    virtual AliEMCALUnfolding *GetClusterUnfoldingInstance() { return fClusterUnfolding; }

  protected:
    AliEMCALGeometry  *fGeom;          // EMCAL geometry
    Float_t            fLogWeight;     // used in AliEMCALRecPoint::EvalGlobalPosition()
    Float_t            fECALocMaxCut;  // this amount of energy must distinguish a local maximum from its neighbours
    Float_t            fMinECut;       // minimum energy of cell   
    TObjArray         *fRecPoints;     //! cluster <=> recPoint
    TClonesArray      *fDigitsArr;     //->   cell <=> digit

    AliEMCALUnfolding *fClusterUnfolding;  // unfolding class instance

    ClassDef(AliEMCALAfterBurnerUF,2)
} ;

#endif // AliEMCALAFTERBURNERUF_H
