/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/// \ingroup graphics
/// \class  AliMpDEVisu
/// \brief GUI for drawing detection element segmentation
///
/// \author Ch. Finck

#ifndef ALI_MP_DE_VISU_H
#define ALI_MP_DE_VISU_H

#include <TGFrame.h>
#include "AliMpPlaneType.h"

class TObject;
class TString;
class TList;
class TArrayI;
class TRootEmbeddedCanvas;
class TGComboBox;
class TGMainFrame;
class TGWindow;
class AliMpVPainter;
class TGNumberEntry;
class TGCheckButton;
class TGTextView;
class AliMpSlat;
class AliMpSector;
class AliMpVSegmentation;
class AliMpDDLStore;
class TGTextEntry;

class AliMpDEVisu : public TGFrame {


public:
    AliMpDEVisu(UInt_t w = 1200, UInt_t h = 600);
    virtual ~AliMpDEVisu();


    void   UpdateComboDE();
    Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);
    void   DrawDE();
    void   NextDE();
    void   DrawManuMotif(Bool_t popup = kFALSE);
    void   DrawQuadrant(Option_t* option, Bool_t popup = kFALSE);
    void   DrawSlat(Option_t* option, Bool_t popup = kFALSE);

    void   ResetManu();
    void   UpdateNameView();
    void   PopUpManuMotif(AliMpSlat* slat);
    void   PopUpManuMotif(AliMpSector* sector);

    void   ClosedPopUpMotif(Int_t id);
    void   InfoDE();
    void   DeletePopUp();
    void   SaveLogMessage();
    void   ClearLogMessage();

    void   HandleMovement(Int_t eventType, Int_t eventX, Int_t eventY, TObject* select);

private:

    const TGWindow*    fkMainWindow; //!< main window
    TGMainFrame*       fMain;        //!< main frame
    TRootEmbeddedCanvas* fEcanvas;   //!< canvas for detection elt

    TGComboBox*    fChamberCombo;    //!< chamber botton 
    TGComboBox*    fDECombo;         //!< DE botton
    TGNumberEntry* fNumberEntry;     //!< manu id button
    TGCheckButton* fPlaneButton;     //!< check button for NB plane, defaultwise B plane
    TGTextView*    fNameDEView;      //!< name of the DE
    TGTextView*    fLogMessage;      //!< log message
    TGTextEntry*   fLogFile;         //!< text entry for log file name
    TList          fTrashList;       //!< list of transient windows to delete

    TArrayI        fDEComboIdx;      //!< array for index vs DE id

    AliMp::PlaneType fCurrentPlane;   //!< current plane type
    Int_t            fCurrentDetElem; //!< current DE
    TString          fCurrentDEName;  //!< current DE name

    const AliMpVSegmentation* fSegmentation; //!< segmentation instance
    AliMpDDLStore*            fDDLStore;     //!< DDL Store

    Int_t            fNumberOfPopUp;   //!< number of manu motif popup window open    

    enum {kChamberCombo, kDECombo, kPlaneType, kDEName, kManuEntries, kLogMessage};


    AliMpDEVisu(const AliMpDEVisu& src);
    AliMpDEVisu& operator=(const AliMpDEVisu& src);

    ClassDef(AliMpDEVisu,1)
};
#endif
