#ifndef ALIITSVBEAMTESTITS04_H
#define ALIITSVBEAMTESTITS04_H

////////////////////////////////////////////////
// ITS geometry class and step manager for the//
//   integrated ITS test beam of Nov. 04      //
////////////////////////////////////////////////
#include "AliITSInitGeometry.h"
#include "AliITS.h"

class TGeoVolume;

class AliITSvBeamTestITS04 : public AliITS {
 
 public:
    AliITSvBeamTestITS04();
    AliITSvBeamTestITS04(const char* name,const char *title);
    virtual ~AliITSvBeamTestITS04();
    virtual void SetNumberOfSPD(Int_t nSPD) {fNspd=nSPD;}
    virtual void SetNumberOfSDD(Int_t nSDD) {fNsdd=nSDD;}
    virtual void SetNumberOfSSD(Int_t nSSD) {fNssd=nSSD;}
    Int_t GetNSPD() const {return fNspd;}
    Int_t GetNSDD() const {return fNsdd;}
    Int_t GetNSSD() const {return fNssd;}
    Int_t GetNumberOfSubDet(const TString& det) const;
    virtual void CreateGeometry();
    virtual void CreateMaterials();
    virtual Int_t IsVersion()const{return (Int_t)kvITS04;}//Return Major Version
    virtual Int_t MinorVersion()const{return fMinorVersion;}//Return Minor Version
    //virtual void InitAliITSgeom();
    virtual void Init();
    //virtual void SetDefaults();
    virtual void StepManager();
    //for writing out geometry set .det write 
    virtual void   SetWriteDet(Bool_t det=kTRUE){fGeomDetOut = det;}
    virtual void   SetWriteDet(const char *f){fWrite[59]=0;strncpy(fWrite,f,59);
                             fGeomDetOut = kTRUE;}// set write file 
    //for reading geometry (JC) set .det read
    virtual void   SetReadDet(Bool_t det=kTRUE){fGeomDetIn = det;}
    virtual void   SetReadDet(const char *f){fRead[59]=0;strncpy(fRead,f,59);
                         fGeomDetIn = kTRUE;} // set read file
 public:
    void AddSPDGeometry(TGeoVolume *moth) const;
    void AddSDDGeometry(TGeoVolume *moth) const;
    void AddSSDGeometry(TGeoVolume *moth) const;
 protected:
           //! ITS mother volume
    static const Int_t fgkNumberOfSPD;       //number of SPD
    static const Int_t fgkNumberOfSDD;       //number of SDD
    static const Int_t fgkNumberOfSSD;       //number of SSD
    // SPD specific
    static const char*    fgSPDsensitiveVolName;  //SPD volume name
    static const Double_t fgkSPDthickness;        //SPD thickness
    static const Double_t fgkSPDwidth;            //SPD width
    static const Double_t fgkSPDlength;           //SPD length
    static const Double_t fgkSPDthickSens;        //SPD sensitive thickness
    static const Double_t fgkSPDwidthSens;        //SPD sensitive width
    static const Double_t fgkSPDlengthSens;       //SPD sensitive length
    static const Double_t fgkSPD0y;               //SPD position
    static const Double_t fgkSPD1y;               //SPD position
    // SDD specific
    static const char*    fgSDDsensitiveVolName;  //SDD volume name
    static const Double_t fgkSDDthickness;        //SDD thickness
    static const Double_t fgkSDDwidth;            //SDD width
    static const Double_t fgkSDDlength;           //SDD length
    static const Double_t fgkSDDthickSens;        //SDD sensitive thickness
    static const Double_t fgkSDDwidthSens;        //SDD sensitive width
    static const Double_t fgkSDDlengthSens;       //SDD sensitive length
    static const Double_t fgkSDD0y;               //SDD position
    static const Double_t fgkSDD1y;               //SDD position
    // SDD specific
    static const char*    fgSSDsensitiveVolName;   //SSD volume name
    static const Double_t fgkSSDthickness;         //SSD thickness
    static const Double_t fgkSSDwidth;             //SSD width
    static const Double_t fgkSSDlength;            //SSD length
    static const Double_t fgkSSDthickSens;         //SSD sensitive thickness
    static const Double_t fgkSSDwidthSens;         //SSD sensitive width
    static const Double_t fgkSSDlengthSens;        //SSD sensitive length
    static const Double_t fgkSSD0y;                //SSD position
    static const Double_t fgkSSD1y;                //SSD position
    //
 private:
    //Int_t GetCurrentLayLaddDet(Int_t &lay,Int_t &ladd, Int_t &det) const;:
    AliITSvBeamTestITS04(const AliITSvBeamTestITS04 &source);//Copy constructor
    //=operator
    AliITSvBeamTestITS04& operator=(const AliITSvBeamTestITS04 &source);
    //
    TGeoVolume *fITSmotherVolume;       //Pointer to mother volume 
    Int_t     fNspd;                    //Number of SPD modules
    Int_t     fNsdd;                    //Number of SDD modules
    Int_t     fNssd;                    //Number of SSD modules
    //for writing out geometry
    Bool_t fGeomDetOut;       // Flag to write .det file out 
    Bool_t fGeomDetIn;        // Flag to read geometry file (JC)
    char   fWrite[60];        //! file name to write .det file 
    char   fRead[60];         // file name to read .det file (JC)
    AliITSVersion_t fMajorVersion; // Major Version
    Int_t           fMinorVersion; // Minor Version
    AliITSInitGeometry fIgm; //! Init geometry object

    ClassDef(AliITSvBeamTestITS04,3) // Full ITS test beam geometry 2004
};

#endif

    
