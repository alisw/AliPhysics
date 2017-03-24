#ifndef UtilMessages_h
#define UtilMessages_h

#include "TError.h"
//Output macros 
const char redTxt[] =     { 0x1b, '[', '1', ';', '3', '1', 'm', 0 };
const char greenTxt[] =   { 0x1b, '[', '1', ';', '3', '2', 'm', 0 };
const char yellowTxt[] =  { 0x1b, '[', '1', ';', '3', '3', 'm', 0 };
const char blueTxt[] =    { 0x1b, '[', '1', ';', '3', '4', 'm', 0 };
const char magentaTxt[] = { 0x1b, '[', '1', ';', '3', '5', 'm', 0 };
const char cyanTxt[] =    { 0x1b, '[', '1', ';', '3', '6', 'm', 0 };
const char normalTxt[] =  { 0x1b, '[', '0', ';', '3', '9', 'm', 0 };

#define Errormsg(where, msg) Error(where, "%s%s%s", redTxt, msg, normalTxt)
#define Warningmsg(where, msg) Warning(where, "%s%s%s", magentaTxt, msg, normalTxt)
#define Infomsg(where, msg) Info(where, "%s%s%s", blueTxt, msg, normalTxt)
#define Infomsgcolor(where, msg, color) Info(where, "%s%s%s", color, msg, normalTxt)
#define Fatalmsg(where, msg) Fatal(where, "%s%s%s", cyanTxt, msg, normalTxt)
#endif
