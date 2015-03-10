/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/// \ingroup mpgraphics
/// \class  AliMpDEVisu
/// \brief GUI for drawing detection element segmentation
///
/// \author Ch. Finck

#ifndef ALI_MP_DE_VISU_H
#define ALI_MP_DE_VISU_H

#include <TGFrame.h>

#include "AliMpPlaneType.h"

#include <TArrayI.h>
#include <TObjArray.h>

class TObject;
class TString;
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
class AliMpManuStore;
class TGTextEntry;
class AliMpMotifPosition;

class AliMpDEVisu : public TGFrame 
{

public:
    AliMpDEVisu(UInt_t w = 1200, UInt_t h = 600);
    virtual ~AliMpDEVisu();


    void   UpdateComboCH();
    void   UpdateComboDE();
    Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);
    void   DrawDE(Bool_t info = kTRUE);
    void   NextDE();
    void   DrawManuMotif(Bool_t popup = kFALSE);
    void   DrawQuadrant(Option_t* option, Bool_t popup = kFALSE);
    void   DrawSlat(Option_t* option, Bool_t popup = kFALSE);

    void   ResetManu();
    void   UpdateNameView(Bool_t firstTime = kFALSE);
    void   PopUpManuMotif(AliMpSlat* slat);
    void   PopUpManuMotif(AliMpSector* sector);
    void   PopUpZoom(Int_t ix0, Int_t iy0, Int_t ix1, Int_t iy1);
    
    void   ClosePopupWindow(Int_t id);
    void   InfoDE();
    void   InfoManuMotif(AliMpMotifPosition* motifPos);
    void   DeletePopUp();
    void   SaveLogMessage();
    void   ClearLogMessage();

    void   HandleMovement(Int_t eventType, Int_t eventX, Int_t eventY, TObject* select);

private:
    void EventToReal(Int_t eventX, Int_t eventY, Double_t& x, Double_t& y) const;
    void CreatePopupWindow(Int_t w, Int_t h, const char* title,
                           AliMpVPainter* painter,
                           const char* option);
    
private:
    /// Not implemented
    AliMpDEVisu(const AliMpDEVisu& src);
    /// Not implemented
    AliMpDEVisu& operator=(const AliMpDEVisu& src);

    const TGWindow*    fkMainWindow; //!<! main window
    TGMainFrame*       fMain;        //!<! main frame
    TRootEmbeddedCanvas* fEcanvas;   //!<! canvas for detection elt

    TGComboBox*    fChamberCombo;    //!<! chamber botton 
    TGComboBox*    fDECombo;         //!<! DE botton
    TGNumberEntry* fNumberEntry;     //!<! manu id button
    TGCheckButton* fPlaneButton;     //!<! check button for NB plane, defaultwise B plane
    TGCheckButton* fZoomButton;      //!<! check button to activate zoom mode, default wise disable
    TGComboBox*    fNameDECombo;     //!<! name of the DE
    TGTextView*    fLogMessage;      //!<! log message
    TGTextEntry*   fLogFile;         //!<! text entry for log file name
    TObjArray      fTrashList;       //!<! list of transient windows to delete

    TArrayI        fDEComboIdx;      //!<! array for index vs DE id
    TString        fNameDEComboIdx[156];  //!<! array for index vs DE names
    TArrayI        fDEOccurrence;     //!<! occurrence of DE

    AliMp::PlaneType fCurrentPlane;   //!<! current plane type
    Int_t            fCurrentDetElem; //!<! current DE
    TString          fCurrentDEName;  //!<! current DE name

    const AliMpVSegmentation* fkSegmentation; //!<! segmentation instance
    AliMpDDLStore*            fDDLStore;      //!<! DDL Store
    AliMpManuStore*           fManuStore;     //!<! Manu Store

    Bool_t           fZoomMode;        //!<! flag for zoom mode on canvas instead of click mode

    enum {kChamberCombo, kDECombo, kPlaneType, kDEName, kManuEntries, kLogMessage, kZoomMode};

    ClassDef(AliMpDEVisu,1) //GUI for drawing detection element segmentation
};
#endif

