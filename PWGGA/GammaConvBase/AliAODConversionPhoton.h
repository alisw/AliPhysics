#ifndef ALIAODCONVERSIONPHOTON_H
#define ALIAODCONVERSIONPHOTON_H

#include "AliConversionPhotonBase.h"
#include "AliKFConversionPhoton.h"
#include "AliAODConversionParticle.h"

class AliAODConversionPhoton : public AliAODConversionParticle, public AliConversionPhotonBase {

  public:

    enum caloPhotonMCFlags_t {
              kIsPhoton               = 0x001, kIsElectron        = 0x002, kIsConversion  = 0x004, kIsConversionFullyContained  = 0x008,
              kIsMerged               = 0x010, kIsMergedPartConv  = 0x020, kIsDalitz      = 0x040, kIsDalitzMerged              = 0x080,
              kIsPhotonWithElecMother = 0x100, kIsShower          = 0x200, kIsSubLeadingEM= 0x400, kIsElectronFromFragPhoton    = 0x800,
    };

    //Constructors
    AliAODConversionPhoton();
    AliAODConversionPhoton(AliKFConversionPhoton *kfphoton);
    AliAODConversionPhoton(TLorentzVector *vec);

    //Copy Constructor
    AliAODConversionPhoton(const AliAODConversionPhoton & g);
    //assignment operator
    AliAODConversionPhoton & operator = (const AliAODConversionPhoton & g);

    //Destructor
    virtual ~AliAODConversionPhoton();

    // Overwrite GetLabelFunctions to Make it accessible via AliAODConversionParticle
    virtual Int_t GetLabel(Int_t i) const { return AliConversionPhotonBase::GetTrackLabel(i); }
    virtual Int_t GetLabel1() const { return AliConversionPhotonBase::GetTrackLabelPositive(); }
    virtual Int_t GetLabel2() const { return AliConversionPhotonBase::GetTrackLabelNegative(); }

    virtual Double_t GetPhotonMass() const {return AliAODConversionParticle::M();}
    virtual Double_t GetPhotonPt() const {return AliAODConversionParticle::Pt();}
    virtual Double_t GetPhotonP() const {return AliAODConversionParticle::P();}
    virtual Double_t GetPhotonEta() const {return AliAODConversionParticle::Eta();}
    virtual Double_t GetPhotonTheta() const {return AliAODConversionParticle::Theta();}
    virtual Double_t GetPhotonPhi() const {return AliAODConversionParticle::Phi();}
    virtual Double_t GetPx() const { return AliAODConversionParticle::Px();}
    virtual Double_t GetPy() const { return AliAODConversionParticle::Py();}
    virtual Double_t GetPz() const { return AliAODConversionParticle::Pz();}
    void CalculateDistanceOfClossetApproachToPrimVtx(const AliVVertex* primVertex);
    void SetMassToZero() { SetE(P()); }

    void SetInvMassPair(Float_t mass) {fInvMassPair=mass;}
    void SetUseForMesonPair(Bool_t useForMeson) {fUseForMesonPair=useForMeson;}
    Float_t GetInvMassPair(){return fInvMassPair;}

    void SetIsTrueConvertedPhoton(){
      fCaloPhoton = 0;
      fCaloPhotonMCFlags = 1;
    }

    Float_t GetDCAzToPrimVtx()const {return fDCAzPrimVtx;}
    Float_t GetDCArToPrimVtx()const {return fDCArPrimVtx;}

    void SetIsCaloPhoton(Int_t clustertype){
        if (clustertype != 4){
            fCaloPhoton = (Char_t)clustertype;
        } else {
            if (GetPhotonPhi() > 1.39 && GetPhotonPhi() < 3.285 )
                fCaloPhoton = 1;
            else 
                fCaloPhoton = 3;
        }
    }
    Char_t GetIsCaloPhoton(){return fCaloPhoton;}
    void SetCaloPhotonMCLabel(Int_t i, Int_t labelCaloPhoton){fCaloPhotonMCLabels[i] = labelCaloPhoton;}
    Int_t GetCaloPhotonMCLabel(Int_t i){return fCaloPhotonMCLabels[i];}
    void SetNCaloPhotonMCLabels(Int_t nLabels){fNCaloPhotonMCLabels = nLabels;}
    Int_t GetNCaloPhotonMCLabels(){return fNCaloPhotonMCLabels;}
    Int_t GetNCaloPhotonMotherMCLabels(){return fNCaloPhotonMotherMCLabels;}
    Int_t GetCaloPhotonMotherMCLabel(Int_t i){return fCaloPhotonMotherMCLabels[i];}
    Int_t GetNNeutralPionMCLabels(){return fNNeutralPionLabels;}
    Int_t GetNeutralPionMCLabel(Int_t i){return fNeutralPionLabels[i];}
    Float_t GetNeutralPionEnergyFraction(Int_t i){return fNeutralPionEnergyFraction[i];}
    Float_t GetLeadingNeutralPionIndex(){return fLeadingNeutralPionIndex;}
    Int_t GetLeadingNeutralPionDaughterIndex(){return fLeadingNeutralPionDaughterIndex;}

    void SetCaloPhotonMCFlags(AliMCEvent *mcEvent, Bool_t enableSort, Bool_t mergedAnalysis = kFALSE, AliVCluster* cluster = NULL);
    void SetCaloPhotonMCFlagsAOD(AliVEvent* event, Bool_t enableSort, Bool_t mergedAnalysis = kFALSE, AliVCluster* cluster = NULL);
    void SetCaloClusterRef(Long_t ref){fCaloClusterRef = ref;}
    Long_t GetCaloClusterRef()const {return fCaloClusterRef;}
    void PrintCaloMCLabelsAndInfo(AliMCEvent *mcEvent);
    void PrintCaloMCFlags ();

    //Calo cluster MC identifiers
    Bool_t IsLargestComponentPhoton(){return fCaloPhotonMCFlags&kIsPhoton;}                       // largest contribution to cluster is photon
    Bool_t IsLargestComponentElectron(){return fCaloPhotonMCFlags&kIsElectron;}                   // largest contribution to cluster is electron
    Bool_t IsConversion(){return fCaloPhotonMCFlags&kIsConversion;}                               // largest contribution to cluster is converted electron
    Bool_t IsConversionFullyContained(){return fCaloPhotonMCFlags&kIsConversionFullyContained;}   // largest contribution to cluster is converted electron & other electron has been found in cluster as well
    Bool_t IsMerged(){return fCaloPhotonMCFlags&kIsMerged;}                                       // largest contribution to cluster is photon, second photon or electron from dalitz decay is found in cluster as well
    Bool_t IsMergedPartConv(){return fCaloPhotonMCFlags&kIsMergedPartConv;}                       // cluster contains more than 1 particle belonging to the same mother particle (i.e. pi0, eta ...)
                                                                                                  // & at least one of the decays was a conversion
    Bool_t IsDalitz(){return fCaloPhotonMCFlags&kIsDalitz;}                                       // cluster contains particle from Dalitz decay
    Bool_t IsDalitzMerged(){return fCaloPhotonMCFlags&kIsDalitzMerged;}                           // cluster contains particle from Dalitz decay & more than one particle of this decay is contained in cluster
    Bool_t IsPhotonWithElecMother(){return fCaloPhotonMCFlags&kIsPhotonWithElecMother;}           // largest contribution to cluster is photon which stems from an electron (i.e. radiation)
    Bool_t IsShower(){return fCaloPhotonMCFlags&kIsShower;}                                       // largest contribution to cluster seems to stem from a shower
    Bool_t IsEMNonLeading(){return !(fCaloPhotonMCFlags&kIsPhoton) && !(fCaloPhotonMCFlags&kIsElectron);} // largest contribution is from hadron
    Bool_t IsSubLeadingEM(){return fCaloPhotonMCFlags&kIsSubLeadingEM;}                           // cluster contains at least one electron or photon from a pi0 or eta in subleading contribution

    Bool_t IsElectronFromFragPhoton(){return fCaloPhotonMCFlags&kIsElectronFromFragPhoton;}      // largest contribution to cluster is from electron, but the converted photon is or stems from fragmentation photon ( q -> q gamma)
    Bool_t IsTrueConvertedPhoton() {
      if (fCaloPhoton == 0 && fCaloPhotonMCFlags == 1) return kTRUE;
        else return kFALSE;
    }
    Bool_t GetUseForMesonPair() {return fUseForMesonPair;}

    Long_t fCaloPhotonMCLabels[50];         //!
    Long_t fCaloPhotonMotherMCLabels[20];   //!
    Long_t fNeutralPionLabels[20];          //!
    Int_t fLeadingNeutralPionIndex;         //!
    Int_t fLeadingNeutralPionDaughterIndex; //!
    Float_t fNeutralPionEnergyFraction[20]; //!
    Long_t fCaloClusterRef;                 //!
    Float_t fDCArPrimVtx;                   //!
    Float_t fDCAzPrimVtx;                   //!
    Float_t fInvMassPair;                   //!
    Int_t fNCaloPhotonMCLabels;             //!
    Int_t fNCaloPhotonMotherMCLabels;       //!
    Int_t fNNeutralPionLabels;              //!
    Int_t fCaloPhotonMCFlags;               //!
    Int_t fPairedId;                        //!
    Char_t fCaloPhoton;                     //!
    Bool_t fUseForMesonPair;                //!

    ClassDef(AliAODConversionPhoton,10)
};


#endif



