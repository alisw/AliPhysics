#ifndef ALIEMCALFOLDER_H
#define ALIEMCALFOLDER_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* History of cvs commits:
 *
* $Log$
* Revision 1.3  2007/10/16 14:36:39  pavlinov
* fixed code violation (almost)
*
* Revision 1.2  2007/09/11 19:38:15  pavlinov
* added pi0 calibration, linearity, shower profile
*
*/

//_________________________________________________________________________
//  Top EMCAL folder - keep everyrhing for calibration task
//  Initial version was created with TDataSet staf
//  TObjectSet -> TFolder; Sep 5, 2007
//                  
//*-- Author: Aleksei Pavlinov (WSU, Detroit, USA) 

// --- ROOT system ---

#include <TFolder.h>
#include <TString.h>

class AliEMCALGeometry;
class AliEMCALSuperModule;
class AliEMCALCell;
class AliESDCaloCluster;
class AliEMCALPi0SelectionParam;
class AliEMCALPi0SelectionParRec;
class AliEMCALCalibData;
class AliEMCALCalibCoefs;
class AliEMCALRecPoint;

class TList;
class TNtuple;

class AliEMCALFolder : public TFolder {

 public:
  
  AliEMCALFolder(); 
  AliEMCALFolder(const AliEMCALFolder& folder); //copy constructor
  AliEMCALFolder(const char* name, const char* title="Top EMCAL folder", Bool_t putToBrowser=kFALSE);
  AliEMCALFolder(const Int_t it, const char* title="Top EMCAL folder", Bool_t putToBrowser=kFALSE);

  virtual ~AliEMCALFolder();

  AliEMCALFolder & operator = (const AliEMCALFolder  & /*rvalue*/) {
    // assignement operator requested by coding convention but not
    // needed                           
    Fatal("operator =", "not implemented");
    return *this;
  };

  void Init(Bool_t putToBrowser=kFALSE);
  // Get methods
  Int_t GetIterationNumber() const {return fCounter;}
  AliEMCALSuperModule* GetSuperModule(const Int_t nm);
  TList*        GetHists() {return fLhists;}
  AliEMCALCell** GetListOfCells() {return fLofCells;}
  AliEMCALCell* GetCell(const Int_t absId);
  void          SetCell(AliEMCALCell *cell, const Int_t absId); 
  AliEMCALPi0SelectionParam*  GetPi0SelectionPar() {return fPi0SelPar;}
  AliEMCALPi0SelectionParRec* GetPi0SelectionParRow(Int_t nrow);
 
  void FillPi0Candidate(const Double_t mgg, AliESDCaloCluster* cl1, AliESDCaloCluster* cl2);
  void FillPi0Candidate(const Double_t mgg, Int_t absIdMax, Int_t nm);
  // Define CC
  void FitAllSMs();    // SM0 now
  // Service routine 
  AliEMCALCalibCoefs* GetCCTable(const char* name) const;
  AliEMCALCalibCoefs* GetCCFirst() {return GetCCTable(fgkCCFirstName.Data());}
  AliEMCALCalibCoefs* GetCCIn() {return GetCCTable(fgkCCinName.Data());}
  AliEMCALCalibCoefs* GetCCOut(){return GetCCTable(fgkCCoutName.Data());}
  Int_t GetSMNumber(AliESDCaloCluster* cl);
  // Recalibration staf - Jun 18,2007
  static AliEMCALRecPoint *GetRecPoint(AliESDCaloCluster *cl,AliEMCALCalibCoefs *tOld,AliEMCALCalibCoefs *tNew, 
  TList *l=0, Double_t deff=-1., Double_t w0=-1., Double_t phiSlope=0.0); 
  // MENU
  //void   Save(const char *fn = "EMCALFOLDER.root", const char *opt="RECREATE");  // *MENU*
  static AliEMCALFolder*  ReadFolder(const char *fn = "EMCALFOLDER.root", const char *opt="READ");
  void   InitAfterRead();                // *MENU*
  void   DrawQA(const int nsm=0); // *MENU*
  void   CreateCellNtuple();      // *MENU*
  void   CreateAndFillAdditionalHists(); // *MENU*

  static const TString GetBaseFolderName() {return fgkBaseFolderName;}
  static const TString GetCCinName()       {return fgkCCinName;}
 protected:
  TList* BookHists();
  Int_t fCounter; // Counter of iteration 
 //
  AliEMCALGeometry *fGeometry; //! pointer to EMCAL geometry
  //
  Int_t     fNumOfCell;     // number of cells as in geometry

  TList*    fLhists;        //! for speed 
  AliEMCALCell** fLofCells; //! unifrom array of cells for fast access; invisible from browser
  AliEMCALPi0SelectionParam* fPi0SelPar; // pi0 selection parameters 
  AliEMCALCalibData *fCalibData; //!
  // 
  TNtuple *fCellNtuple;  //! for quick cell anaylsis

  static const TString fgkBaseFolderName;  // base name of EMCAL Folder  
  static const TString fgkCCFirstName;     // name of first calib.table 
  static const TString fgkCCinName;        // name of initial calib.coefs. table 
  static const TString fgkCCoutName;       // name of out calib.coefs. table 
  static const TString fgkDirOfRootFiles;  // name of directory for saving EMCAL folder

  TList *fLobj; // list of all objects

  void TestSMStruct();  // *MENU*

  ClassDef(AliEMCALFolder,3) // EMCAL folder
    
};

#endif // ALIEMCALFOLDER_H
