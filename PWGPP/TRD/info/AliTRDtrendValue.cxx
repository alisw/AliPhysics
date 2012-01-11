////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Trend Value Incapsulation                                             //
//                                                                        //
//  Authors:                                                              //
//    Alexandru Bercuci <A.Bercuci@gsi.de>                                //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"

#include "AliLog.h"
#include "AliTRDtrendValue.h"

ClassImp(AliTRDtrendValue)


//____________________________________________
AliTRDtrendValue::AliTRDtrendValue() 
  : TNamed("none", "none")
  ,fAlarmLevel(0)
  ,fValue(0.)
  ,fResponsible()
  ,fNnotifiable(0)
{
//  Constructor. Reset all fields.
  memset(fLimits, 0, 2*(kNlevels+1)*sizeof(Double_t));
  for(Int_t ilevel(kNlevels); ilevel--; ) snprintf(fAlarmMessage[ilevel], 1024, " ");
}

//____________________________________________
AliTRDtrendValue::AliTRDtrendValue(Char_t *n, Char_t *t) 
  : TNamed("none", t)
  ,fAlarmLevel(0)
  ,fValue(0.)
  ,fResponsible()
  ,fNnotifiable(0)
{
//  Constructor. Define name and title for trend variable.
  TString s(n);
  TObjArray *names(s.Tokenize("_"));
  if(names->GetEntriesFast()!=2){
    AliError(Form("Wrong trend value name format. Trend value name should be of the form \"trendClass_trendValue\" with only one \"_\" character."));
  } else SetName(n);

  memset(fLimits, 0, 2*(kNlevels+1)*sizeof(Double_t));
  for(Int_t ilevel(kNlevels); ilevel--; ) snprintf(fAlarmMessage[ilevel], 1024, " ");
}

//____________________________________________
Int_t AliTRDtrendValue::GetAlarmLevel()
{
  // check value against limits and do some more work
  fAlarmLevel=kNlevels-1;
  for(Int_t ilevel(0); ilevel<kNlevels+1; ilevel++)
    if(fValue<fLimits[2*ilevel+1] &&
       fValue>=fLimits[2*ilevel]){ 
      fAlarmLevel = ilevel;
      break;
    }

  return fAlarmLevel;
}

//____________________________________________
const char* AliTRDtrendValue::GetAlarmMessage() const
{
// Check if value triggered alarm
  if(!fAlarmLevel) return "OK";
  else return fAlarmMessage[fAlarmLevel-1];
}

//____________________________________________
const char* AliTRDtrendValue::GetClassName() const
{
// Check task to which value belong
  TString s(TNamed::GetName());
  TObjArray *names(s.Tokenize("_"));
  if(names->GetEntriesFast()!=2){
    AliError(Form("Wrong trend value name format."));
    return NULL;
  }

  return ((TObjString*)names->At(0))->String().Data();
}

//____________________________________________
const char* AliTRDtrendValue::GetValueName() const
{
// value name
  TString s(TNamed::GetName());
  TObjArray *names(s.Tokenize("_"));
  if(names->GetEntriesFast()!=2){
    AliError(Form("Wrong trend value name format."));
    return NULL;
  }
  return ((TObjString*)names->At(1))->String().Data();
}

//____________________________________________
const char* AliTRDtrendValue::GetResponsible(Char_t *n, Char_t *mail) const
{
// Get responsible with name and mail
  if(n) snprintf(n, 100, "%s", fResponsible.fNameR);
  if(mail) snprintf(mail, 200, "%s", fResponsible.fMail);
  return Form("%s <%s>", fResponsible.fNameR, fResponsible.fMail);
}

//____________________________________________
const char* AliTRDtrendValue::GetNotifiable(Int_t in, Char_t *n, Char_t *mail) const
{
// Get noticible person "in" with name and mail
  if(in<0||in>=fNnotifiable) return NULL;
  if(n) snprintf(n, 100, "%s", fNotifiable[in].fNameR);
  if(mail) snprintf(mail, 200, "%s", fNotifiable[in].fMail);
  return Form("%s <%s>", fNotifiable[in].fNameR, fNotifiable[in].fMail);
}

//____________________________________________
void AliTRDtrendValue::SetNotifiable(const Char_t *name, const Char_t *mail)
{
// add noticible person to DB
  if(fNnotifiable==kNnotifiable){
    AliWarning(Form("Could not add %s for notification. Only %d persons can be registered for notification.", name, kNnotifiable));
    return;
  }
  snprintf(fNotifiable[fNnotifiable].fNameR, 100, "%s", name);
  snprintf(fNotifiable[fNnotifiable].fMail, 200, "%s", mail);
  fNnotifiable++;
}

//____________________________________________
void AliTRDtrendValue::SetResponsible(const Char_t *name, const Char_t *mail) 
{
// set responsible person for trend
  snprintf(fResponsible.fNameR, 100, "%s", name);
  snprintf(fResponsible.fMail, 200, "%s", mail);
}

//____________________________________________
void AliTRDtrendValue::Print(Option_t */*o*/) const
{
//   name - title
//   value - limits
//   alarm level, message
//   responsible

  printf("    %s [%s] - %s\n", GetValueName(), GetClassName(), GetTitle());
  printf("*** %f limits[%f %f]\n", fValue, fLimits[0], fLimits[1]);
  if(fAlarmLevel){
    printf("*** Alarm level   : %d limits[%f %f]\n", fAlarmLevel, fLimits[2*fAlarmLevel], fLimits[2*fAlarmLevel+1]);
    printf("*** Alarm message : %s\n", GetAlarmMessage());
  }
  printf("*** Responsible %s <%s>\n", fResponsible.fNameR, fResponsible.fMail);
  if(fNnotifiable){
    printf("*** Notifiable person(s) ***\n");
    for(Int_t i(0); i<fNnotifiable; i++)
      printf("        %s <%s>\n", fNotifiable[i].fNameR, fNotifiable[i].fMail);
  }
}

//____________________________________________
AliTRDtrendValue::AliTRDtrendValueResponsible::AliTRDtrendValueResponsible(Char_t *n, Char_t *m) 
{
// define person with mail and mail
  if(n) snprintf(fNameR, 100, "%s", n); else snprintf(fNameR, 100, " ");
  if(m) snprintf(fMail, 200, "%s", m); else snprintf(fMail, 200, " ");
}
