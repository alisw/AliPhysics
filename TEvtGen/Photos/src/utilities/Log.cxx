#include <fstream>
#include "Log.h"
using std::streambuf;
using std::stringstream;
using std::ostream;
using std::cout;
using std::cerr;
using std::endl;

namespace Photospp
{

void (*PHOERR)(int,const char*,double) = Log::PHOERR;
void (*PHOREP)()                       = Log::PHOREP;

list<Log::Pointer*> *Log::PointerList = NULL;

streambuf   *Log::bCout=cout.rdbuf(),*Log::bCerr=cerr.rdbuf();
ostream     *Log::out=&cout;
stringstream Log::buf;
int  Log::warnLimit=100;
int  Log::decays[4] = {0};
int  Log::dCount =0,Log::dRangeS =65535,Log::dRangeE =65534;
int  Log::faCount=0,Log::faRangeS=65535,Log::faRangeE=65534;
int  Log::iCount =0,Log::wCount =0,Log::eCount =0,Log::asCount=0, Log::asFailedCount=0;
bool Log::iAction=1,Log::wAction=1,Log::eAction=1,Log::asAction=1,Log::rAction=1;

void Log::AddDecay(int type)
{
	decays[type]++;
}

ostream& Log::Debug(unsigned short int code, bool count)
{
	if(count) ++dCount;
	if(code>=dRangeS && code<=dRangeE ) return *out<<"DEBUG("<<code<<") from PHOTOS:"<<endl;
	return buf.seekp(0);
}


ostream& Log::Info(bool count)
{
	if(count) ++iCount;
	if(iAction) return *out<<"INFO from PHOTOS:"<<endl;
	return buf.seekp(0);
}


ostream& Log::Warning(bool count)
{
	if(count) ++wCount;
	if(warnLimit>0 && wCount>=warnLimit)
	{
		if(wAction)
		{
			*out<<"WARNING from PHOTOS:"<<endl<<"Limit reached ("<<warnLimit<<"). Warnings suppressed."<<endl;
			wAction=false;
		}
		return buf.seekp(0);
	}
	if(wAction && count) return *out<<"WARNING from PHOTOS:"<<endl;
	if(wAction)          return *out;
	return buf.seekp(0);
}


ostream& Log::Error(bool count)
{
	if(count) ++eCount;
	if(eAction) return *out<<"ERROR from PHOTOS:"<<endl;
	buf.seekp(0);
	return buf;
}

void Log::Assert(bool check, char *text)
{
	++asCount;
	if(check) return;
	++asFailedCount;
	if(text==NULL)	*out<<"ASSERT from PHOTOS:"<<endl<<"Assertion failed. "<<endl;
	else *out<<"ASSERT from PHOTOS:"<<endl<<"Assertion failed: "<<text<<endl;
	if(asAction) exit(-1);
}

void Log::Fatal(string text,unsigned short code)
{
	++faCount;
	if(text.size()==0) *out<<"FATAL ERROR from PHOTOS:"<<endl<<"Terminated by a call to Log::Exit();"<<endl;
	else *out<<"FATAL ERROR from PHOTOS: "<<endl<<text<<endl;
	if(code<faRangeS || code>faRangeE) exit(-1);
}

void Log::RedirectOutput(void (*func)(), ostream& where)
{

	if(!rAction) { func(); return; }
	cout.rdbuf(where.rdbuf());
	cerr.rdbuf(where.rdbuf());
	where<<endl;
	func();
	cout.rdbuf(bCout);
	cerr.rdbuf(bCerr);
}

void Log::RedirectOutput(ostream& where)
{
	if(!rAction) return;
	cout.rdbuf(where.rdbuf());
	cerr.rdbuf(where.rdbuf());
	where<<endl;
}

void Log::Summary()
{
	*out<<"---------------------------- Photos Log Summary ------------------------------"<<endl;
	*out<<" Debug:   \t";
	if(dRangeS>dRangeE) *out<<"(OFF)";
	*out<<"\t\t"<<dCount<<"\t";
	if(dRangeS<=dRangeE) *out<<"Debug range: "<<dRangeS<<" - "<<dRangeE;
	*out<<endl;
	*out<<" Info:    \t";
	if(!iAction) *out<<"(OFF)";
	*out<<"\t\t"<<iCount<<"\t"<<endl;
	*out<<" Warnings:\t";
	if(!wAction) if(warnLimit>0 && wCount>warnLimit) *out<<"(SUPP.)"; else *out<<"(OFF)";
	*out<<"\t\t"<<wCount<<"\t"<<endl;
	*out<<" Errors:  \t";
	if(!eAction) *out<<"(OFF)";
	*out<<"\t\t"<<eCount<<"\t"<<endl;
	if(asCount || !asAction || faRangeS<faRangeE) cout<<"-----------------------------------"<<endl;
	if(asCount>0) *out<<" Asserts:\t\t\t"<<asCount<<endl;
	if(!asAction) *out<<" Failed asserts ignored:\t"<<asFailedCount<<endl;
	if(faRangeS<=faRangeE) *out<<" Fatal errors ignored:  \t"<<faCount<<endl;
	cout<<"-----------------------------------"<<endl;
	if(decays[3]) cout<<" Normal decays:                        "<<decays[3]<<endl;
	if(decays[2]) cout<<" Decays without mother:                "<<decays[2]<<endl;
	if(decays[1]) cout<<" Decays without mother & grandmothers: "<<decays[1]<<endl;
	if(decays[0]) cout<<" Decayed using Tauola gun:             "<<decays[0]<<endl;
	*out<<"------------------------------------------------------------------------------"<<endl;
}


//----------------------------------------------------------------------
//
//    PHOTOS:   PHOton radiation in decays ERRror handling
//
//    Purpose:  Inform user  about (fatal) errors and warnings generated
//              by either the user or the program.
//
//    Input Parameters:   IMES, TEXT, DATA
//
//    Output Parameters:  None
//
//    Author(s):  B. van Eijk                     Created at:  29/11/89
//                                                Last Update: 18/06/13
//
//----------------------------------------------------------------------
void Log::PHOERR(int IMES,const char *TEXT,double DATA){

  static int IERROR=0;
  double  SDATA;
  static int PHOMES=10;
  static int i=1;
  char star80[81]= "********************************************************************************";

  if (IMES<=PHOMES) phosta_.status[IMES-i]=phosta_.status[IMES-i]+1;
// 
//    Count number of non-fatal errors...
  if ((IMES ==  6) && (phosta_.status[IMES-i]>=2)) return;
  if ((IMES == 10) && (phosta_.status[IMES-i]>=2)) return;
  SDATA=DATA;
  //  int PHLUN=(int)pholun_.phlun;
  bool IFSTOP=phosta_.ifstop;
  FILE *PHLUN = stdout;
  int furthA=0;
  fprintf(PHLUN,"%s\n",star80);
  fprintf(PHLUN,"*\n");  //9120
  //      GOTO (10,20,30,40,50,60,70,80,90,100),IMES

  switch(IMES){
  case 1:
    fprintf(PHLUN,"* %s: Too many charged Particles, NCHARG = %6i\n", TEXT,(int)SDATA);   //I6
    furthA= 110;
    break;
  case 2:
    fprintf(PHLUN,"* %s: Too much Bremsstrahlung required, PRSOFT = %15.6f\n", TEXT,SDATA);//F15.6
    furthA= 110;
    break;
  case 3:
    fprintf(PHLUN,"* %s: Combined Weight is exceeding 1., Weight = %15.6f\n", TEXT,SDATA);   //F15.6
    furthA= 110;
    break;
  case 4:
    fprintf(PHLUN,"* %s: Error in Rescaling charged and neutral Vectors\n", TEXT);
    furthA= 110;
    break;
  case 5:
    fprintf(PHLUN,"* %s: Non matching charged Particle Pointer, NCHARG = %5i\n", TEXT,(int)SDATA);  //I5
    furthA= 110;
    break;
  case 6:
    fprintf(PHLUN,"* %s: Do you really work with a Particle of Spin: %4.1f\n", TEXT,SDATA);   //F4.1
    furthA= 130;
    break;
  case 7:
    fprintf(PHLUN,"* %s: Stack Length exceeded, NSTACK = %5i\n", TEXT,(int)(SDATA));//I5
    furthA= 110;
    break;
  case 8:
    fprintf(PHLUN,"* %s: Random Number Generator Seed(1) out of Range: %8i\n", TEXT,(int)SDATA);//I8
    furthA= 110;
    break;
  case 9:
    fprintf(PHLUN,"* %s: Random Number Generator Seed(2) out of Range: %8i\n", TEXT,(int)SDATA);//I8
    furthA= 110;
    break;
  case 10:
    fprintf(PHLUN,"* %s: Available Phase Space below Cut-off: %15.6f GeV/c^2\n", TEXT,SDATA);//F15.6
    furthA= 130;
    break;
  default:
    fprintf(PHLUN,"* Funny Error Message: %4i ! What to do ?\n", IMES);//I4
    furthA= 120;
    break;
  }

 switch(furthA){
 case 110:
   fprintf(PHLUN,"* Fatal Error Message, I stop this Run !\n");
   fprintf(PHLUN,"*\n"); //9120
   fprintf(PHLUN,"%s\n",star80);
   if (IFSTOP){ 
     exit(-1);
   }
   else{
     fprintf(PHLUN,"*\n"); //9120
     fprintf(PHLUN,"%s\n",star80);
     break;
   }      
 case 120:
   IERROR=IERROR+1;
   if (IERROR>=10){
     fprintf(PHLUN,"* 10 Error Messages generated, I stop this Run !\n");
     fprintf(PHLUN,"*\n");//9120
     fprintf(PHLUN,"%s\n",star80);
     if (IFSTOP){
       exit(-1);
     }
     else{
       fprintf(PHLUN,"*\n"); //9120
       fprintf(PHLUN,"%s\n",star80);
       break;
     }
   }  
 case 130:
  fprintf(PHLUN,"*\n");  //9120
  fprintf(PHLUN,"%s\n",star80);
  break;
 }
 return;


 //9120 FORMAT(1H ,'*',T81,'*')
 // 9140 FORMAT(1H ,'* Fatal Error Message, I stop this Run !',T81,'*')
 // 9150 FORMAT(1H ,'* 10 Error Messages generated, I stop this Run !',T81,
 //     &'*')
}


//----------------------------------------------------------------------
//
//    PHOTOS:   PHOton radiation in decays run summary REPort
//
//    Purpose:  Inform user about success and/or restrictions of PHOTOS
//              encountered during execution.
//
//    Input Parameters:   Common /PHOSTA/
//
//    Output Parameters:  None
//
//    Author(s):  B. van Eijk                     Created at:  10/01/92
//                                                Last Update: 18/06/13
//
//----------------------------------------------------------------------
void Log::PHOREP(){
  static int PHOMES=10;
  int I;
  bool ERROR=false;
  //  int PHLUN=(int)pholun_.phlun;
  char star80[81]= "********************************************************************************";
  char X26[27] = "                          ";
  char EQ25[26]= "=========================";
  char X30[31] = "                              ";
  char X22[23] = "                      ";
  char X23[24 ]= "                       ";
  char X16[17] = "                ";
  FILE *PHLUN = stdout;
  fprintf(PHLUN," \n");
  fprintf(PHLUN,"%s\n",star80);
  fprintf(PHLUN,"*\n");
  fprintf(PHLUN,"* %s %s\n",X26,EQ25);
  fprintf(PHLUN,"* %s PHOTOS Run Summary\n",X30);
  fprintf(PHLUN,"* %s %s\n",X26,EQ25);
  fprintf(PHLUN,"*\n");
  for(I=1;I<=PHOMES;I++){

    if (phosta_.status[I-1] == 0) break;
    if ((I == 6)|| (I == 10)){
      fprintf(PHLUN,"* %s Warning # %2i  occured %6i times\n",X22, I,phosta_.status[I-1]); // I2 I6 
    }
    else{
      ERROR=true;
      fprintf(PHLUN,"* %s Error # %2i occured %6i  times\n",X23, I,phosta_.status[I-1]);// I2 I6
    }	      
  }

  if (!ERROR) fprintf(PHLUN,"* %s PHOTOS Execution has successfully terminated\n",X16);
  fprintf(PHLUN,"*\n");
  fprintf(PHLUN,"%s\n",star80);
  return;

//      RETURN
// 9000 FORMAT(1H1)
// 9010 FORMAT(1H ,80('*'))
// 9020 FORMAT(1H ,'*',T81,'*')
// 9030 FORMAT(1H ,'*',26X,25('='),T81,'*')
// 9040 FORMAT(1H ,'*',30X,'PHOTOS Run Summary',T81,'*')
// 9050 FORMAT(1H ,'*',22X,'Warning #',I2,' occured',I6,' times',T81,'*')
// 9060 FORMAT(1H ,'*',23X,'Error #',I2,' occured',I6,' times',T81,'*')
// 9070 FORMAT(1H ,'*',16X,'PHOTOS Execution has successfully terminated',
//     &T81,'*')
}





} // namespace Photospp
