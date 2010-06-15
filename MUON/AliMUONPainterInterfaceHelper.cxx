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

// $Id$

#include "AliMUONPainterInterfaceHelper.h"

///\class AliMUONPainterInterfaceHelper
///
/// Helper class to work with TGButtonGroup
///
/// This class only works if the buttons in the TGButtonGroup have contiguous
/// ids, and if those ids start at ButtonStartingId().
/// Not bullet-proof, I admit, but this is the only way I've found with
/// the current TGButtonGroup implementation which does not allow to loop
/// easily on all buttons.
///
// Author Laurent Aphecetche, Subatech

#include "AliMUONPainterEnv.h"
#include "AliMUONPainterHelper.h"
#include "AliLog.h"
#include <Riostream.h>
#include <TClass.h>
#include <TGButtonGroup.h>
#include <TGButton.h>
#include <TGClient.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <RVersion.h>
#include <cassert>

///\cond CLASSIMP
ClassImp(AliMUONPainterInterfaceHelper)
///\endcond

//_____________________________________________________________________________
AliMUONPainterInterfaceHelper::AliMUONPainterInterfaceHelper() : TObject()
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONPainterInterfaceHelper::~AliMUONPainterInterfaceHelper()
{
  /// dtor
}

//_____________________________________________________________________________
void
AliMUONPainterInterfaceHelper::AddRadioButton(TGButtonGroup& bg,
                                              const TString& name,
                                              void* userData,
                                              Bool_t select)
{
  /// Add a radio button to a group
  Int_t n = bg.GetCount();
  
  TGButton* button = new TGRadioButton(&bg,
                                       name.Data(),
                                       n+ButtonStartingId());
  button->SetUserData(userData);
  button->SetOn(select,kFALSE);  
}

//_____________________________________________________________________________
void
AliMUONPainterInterfaceHelper::AddCheckButton(TGButtonGroup& bg,
                                              const TString& name,
                                              void* userData,
                                              Bool_t select)
{
  /// Add a check button to a group
  Int_t n = bg.GetCount();
  
  TGButton* button = new TGCheckButton(&bg,
                                       name.Data(),
                                       n+ButtonStartingId());
  button->SetUserData(userData);
  button->SetOn(select,kFALSE);
}

//_____________________________________________________________________________
void 
AliMUONPainterInterfaceHelper::ClearButtons(TGButtonGroup& bg)
{
  /// Remove all buttons from group.
  /// Bear in mind that you're thus disconnecting the group from
  /// any signals it might have. So you must re-connect them after
  /// a call to this method, if you want to.
  
  while ( bg.GetCount() > 0 )
  {
    TGTextButton* button = (TGTextButton*)(bg.GetButton(bg.GetCount()));
    if ( !button ) 
    {
      AliFatalClass(Form("Got a null button for bg.Count()=%d",bg.GetCount()));
    }
    bg.Remove(button);
#if ROOT_VERSION_CODE <= ROOT_VERSION(5,16,0)
    button->DestroyWindow();
#endif
    delete button;
  }
}

//_____________________________________________________________________________
void
AliMUONPainterInterfaceHelper::Copy(const TGButtonGroup& src, TGButtonGroup& dest)
{
  /// Copy a button group into another one
  AliDebugClass(1,Form("src=%p (%s) count=%d dest=%p (%s) count=%d",
                       &src,src.GetTitle(),src.GetCount(),
                       &dest,dest.GetTitle(),dest.GetCount()));

  StdoutToAliDebugClass(1,cout << "---copying:" << endl; Dump(src);
                        cout << "---to:" << endl; Dump(dest));
  
  ClearButtons(dest);
  
  dest.SetTitle(src.GetTitle());
  
  if ( &src != &dest )
  {
    for ( Int_t i = ButtonStartingId(); i < ButtonStartingId() + src.GetCount(); ++i )
    {
      TGTextButton* tb = static_cast<TGTextButton*>(src.GetButton(i));
      TGButton* button = new TGRadioButton(&dest,tb->GetTitle(),tb->WidgetId());
      assert(tb->WidgetId()==i);
      button->SetUserData(tb->GetUserData());
      button->SetOn(tb->IsOn(),kFALSE);
    }
  }
}

//_____________________________________________________________________________
void 
AliMUONPainterInterfaceHelper::Dump(const TGButtonGroup& bg)
{
  /// Printout
  cout << Form("ButtonGroup %s %s",bg.GetName(),bg.GetTitle()) << endl;
  
  for ( Int_t i = ButtonStartingId(); i < bg.GetCount() + ButtonStartingId(); ++i ) 
  {
    TGTextButton* button = static_cast<TGTextButton*>(bg.GetButton(i));
    if ( button ) 
    {
      cout << Form("i %2d button %s id %d wid %d ON %d",
                 i,button->GetTitle(),button->GetId(),
                 button->WidgetId(),
                 button->IsOn()) << endl;
    }
    else
    {
      cout << Form("i %2d button = 0x0",i) << endl;
    }
  }
}

//_____________________________________________________________________________
TGButton* 
AliMUONPainterInterfaceHelper::FindButtonByName(const TGButtonGroup& bg, 
                                                const TString& name)
{
  /// Find a button by name
  
  for ( Int_t i = ButtonStartingId(); i < ButtonStartingId() + bg.GetCount(); ++i )
  {
    TGTextButton* button = static_cast<TGTextButton*>(bg.GetButton(i));
    if (!button)
    {
      AliErrorClass(Form("(%s) : Something wrong",name.Data()));
      Dump(bg);
    }
    else
    {
      if ( name == button->GetTitle() )
      {
        return button;
      }
    }
  }
  return 0x0;
}

//_____________________________________________________________________________
TGButton* 
AliMUONPainterInterfaceHelper::FindButtonByUserData(const TGButtonGroup& bg, 
                                                    const void* userData)
{
  /// Find a button by userData
  
  for ( Int_t i = ButtonStartingId(); i < ButtonStartingId() + bg.GetCount(); ++i )
  {
    TGButton* button = bg.GetButton(i);
    if ( button->GetUserData() == userData )
    {
      return button;
    }
  }
  return 0x0;
}

//_____________________________________________________________________________
TGButton* 
AliMUONPainterInterfaceHelper::FindDownButton(const TGButtonGroup& bg)
{
  /// Find which button is down (in a radio group)
  
  AliDebugClass(1,Form("bg %s",bg.GetTitle()));
  
  for ( Int_t i = ButtonStartingId(); i < ButtonStartingId() + bg.GetCount(); ++i )
  {
    TGButton* button = bg.GetButton(i);
    if ( button->IsOn() ) 
    {
      AliDebugClass(1,Form("button %s",button->GetTitle()));
      return button;
    }
  }
  return 0x0;
}


//_____________________________________________________________________________
void 
AliMUONPainterInterfaceHelper::SetBackgroundColor(const char* resourceBaseName, 
                                                  TGWindow& window)
{
  /// Set the background color of the window
  TString rs(Form("%s.BackgroundColor",resourceBaseName));
  TString colorName = AliMUONPainterHelper::Instance()->Env()->String(rs.Data(),"#c0c0c0");

  Pixel_t color;
  Bool_t ok = gClient->GetColorByName(colorName, color);
  if ( ok ) 
  {
    window.SetBackgroundColor(color);
    AliDebugClass(1,Form("Setting %s color to %s",rs.Data(),colorName.Data()));
  }
}

//_____________________________________________________________________________
void 
AliMUONPainterInterfaceHelper::SetState(TGButtonGroup& bg, Bool_t state)
{
  /// should not be needed with root > 5.16/00 as there's a TGButtonGroup::SetState
  /// method now...

#if ROOT_VERSION_CODE > ROOT_VERSION(5,16,0)
  bg.SetState(state);
#else
  for ( Int_t i = ButtonStartingId(); i < ButtonStartingId() + bg.GetCount(); ++i ) 
  {
    TGTextButton* button = (TGTextButton*)(bg.GetButton(i));
    if ( state ) 
    {
      button->SetState(kButtonUp);
    }
    else
    {
      button->SetState(kButtonDisabled);
    }
  }
#endif  
}

//_____________________________________________________________________________
void 
AliMUONPainterInterfaceHelper::Select(TGButtonGroup& bg, 
                                      const TString& buttonName,
                                      Bool_t emit)
{
  /// Select which button should be on
  
  AliDebugClass(1,Form("bg %s buttonName %s",bg.GetTitle(),buttonName.Data()));
  
  for ( Int_t i = ButtonStartingId(); i < ButtonStartingId() + bg.GetCount(); ++i ) 
  {
    TGTextButton* button = (TGTextButton*)(bg.GetButton(i));
    if ( buttonName == button->GetTitle() || buttonName == "*" ) 
    {
      if ( emit ) 
      {
        button->Clicked();
      }
      else
      {        
        button->SetOn(kTRUE);
      }
    }
  }  
}

//_____________________________________________________________________________
void 
AliMUONPainterInterfaceHelper::Unselect(TGButtonGroup& bg, 
                                        const TString& buttonName,
                                        Bool_t emit)
{
  /// Unselect a given button
  
  for ( Int_t i = ButtonStartingId(); i < ButtonStartingId() + bg.GetCount(); ++i ) 
  {
    TGTextButton* button = (TGTextButton*)(bg.GetButton(i));
    if ( buttonName == button->GetTitle() || buttonName == "*" ) 
    {
      if ( emit ) 
      {
        button->Released();
      }
      else
      {        
        button->SetOn(kFALSE);
      }
    }
  }  
}

//_____________________________________________________________________________
void
AliMUONPainterInterfaceHelper::RemoveButton(TGButtonGroup& bg,
                                            const TGButton* button)
{
  /// Remove a button
  
  // need to redo it from scratch in order not the leave holes in the
  // id numbering
  
  TGButtonGroup bgtmp(bg.GetParent(),bg.GetTitle());
  
  Int_t id(ButtonStartingId());
  
  for ( Int_t i = ButtonStartingId(); i < ButtonStartingId() + bg.GetCount(); ++i ) 
  {
    TGTextButton* b = static_cast<TGTextButton*>(bg.GetButton(i));

    if ( b != button ) 
    {
      TGButton* bc = new TGRadioButton(&bgtmp,b->GetTitle(),id);
      ++id;
      bc->SetUserData(b->GetUserData());
      bc->SetOn(b->IsOn());
    }    
  }    

  ClearButtons(bg);

  Copy(bgtmp,bg);
  
  bg.Show();
}
