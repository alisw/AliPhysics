/******************************************************************************
 *                      T H E R M I N A T O R                                 *
 *                   THERMal heavy-IoN generATOR                              *
 *                           version 1.0                                      *
 *                                                                            *
 * Authors of the model: Wojciech Broniowski, Wojciech.Broniowski@ifj.edu.pl, *
 *                       Wojciech Florkowski, Wojciech.Florkowski@ifj.edu.pl  *
 * Authors of the code:  Adam Kisiel, kisiel@if.pw.edu.pl                     *
 *                       Tomasz Taluc, ttaluc@if.pw.edu.pl                    *
 * Code designers: Adam Kisiel, Tomasz Taluc, Wojciech Broniowski,            *
 *                 Wojciech Florkowski                                        *
 *                                                                            *
 * For the detailed description of the program and furhter references         * 
 * to the description of the model plesase refer to: nucl-th/0504047,         *
 * accessibile at: http://www.arxiv.org/nucl-th/0504047                       *
 *                                                                            *
 * Homepage: http://hirg.if.pw.edu.pl/en/therminator/                         *
 *                                                                            *
 * This code can be freely used and redistributed. However if you decide to   *
 * make modifications to the code, please contact the authors, especially     *
 * if you plan to publish the results obtained with such modified code.       *
 * Any publication of results obtained using this code must include the       *
 * reference to nucl-th/0504047 and the published version of it, when         *
 * available.                                                                 *
 *                                                                            *
 *****************************************************************************/
#include <TMath.h>
#include "THGlobal.h"
#include "ReadPar.h"
#include "Parser.h"
#include "DecayChannel.h"
#include <sstream>
#include <cstring>

extern ReadPar *sRPInstance;

const double factorials[7] = {1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0 };

Parser::Parser(ParticleDB *aDB)
{
  mDB = aDB;
  ReadParameters();
}

Parser::~Parser()
{
}

int Parser::check(char *a,char *b)
{
  int tIter3=0;

  while(a[tIter3]!='\0') 
    {
      if(a[tIter3]!=b[tIter3]) return 0;
      tIter3++;
    }
  return 1;
}

int Parser::GetParticleCount()
{
  return mParticleCount;
}

char * Parser::GetParticleName(int i)
{
  return mNameBuffer[i];
}

double Parser::GetParticleMass(int i)
{
  return mMassBuffer[i];
}

double Parser::GetParticleGamma(int i)
{
  return mGammaBuffer[i];
}

double Parser::GetParticleSpin(int i)
{
  return mSpinBuffer[i];
}

int Parser::GetParticleBarionN(int i)
{
  return mBarionBuffer[i];
}

int Parser::GetParticleStrangeN(int i)
{
  return mStrangeBuffer[i];
}

double Parser::GetParticleI3(int i)
{
  return mI3Buffer[i];
}

int Parser::GetParticleDecayChannelCount(int i,int j)
{
  return mDecayBuffer[i][j];
}

int Parser::GetParticleNumber(int i)
{
  return mTypeCountBuffer[i];
}

void Parser::ReadInput()
{
  int j,tPartIter=0,l,tIter2, tIter; //variables
  char str[200];
  char str1[200];
  double spin1,spin2,value;

  //////  START OF "TABLES.M" /////////
  ifstream in("tables.m");
  if(in)
    {
      //START OF HEAD-LINE
      in.ignore(100,'\n');
      //END OF HEAD-LINE

      for(tIter2=0;tIter2<50;tIter2++) str[tIter2]='\0';

      tPartIter=0;

      //START OF DATA
      while(in>>str)
	{
	  if(*str == 'x') break;

	  l=0;
	  //PARTICLE NAME AND MASS
	  if(str[0] == 'm' && str[1] == '[')
	    {
	      for(;;)
		{
		  if(str[l+2] == ']') break;
		  mNameBuffer[tPartIter][l]=str[l+2];	//name
		  l++;
		}
	      //		mNameBuffer[tPartIter][l]='\0';
	      in>>str;	// sign "="
	      in>>value;	// mass
	      mMassBuffer[tPartIter]=value;
	      mTypeCountBuffer[tPartIter]=tPartIter;
	      for(tIter2=0;tIter2<50;tIter2++) str[tIter2]='\0';
	    }
	    
	  if(str[0] == 'G' && str[3] == '[')
	    {
	      in>>str;	// sign "="
	      in>>value;	// gamma
	      mGammaBuffer[tPartIter]=value;
	      for(tIter2=0;tIter2<50;tIter2++) str[tIter2]='\0';
	    }
	    
	  // spin
	  if(str[0] == 'J' && str[1] == '[')
	    {
	      in>>str;	// sign "="
	      in>>str;
	      if(str[1] == '/')
		{
		  *str1=str[0];spin1=atof(str1);
		  *str1=str[2];spin2=atof(str1);
		  mSpinBuffer[tPartIter]=spin1/spin2;
		}
	      if(str[0] == '-')
		{
		  *str1=str[1];spin1=atof(str1);
		  mSpinBuffer[tPartIter]=-spin1;
		}
	      if(str[0]!='-' && str[1]!='/') 
		{
		  *str1=str[0];
		  spin1=atof(str1);
		  mSpinBuffer[tPartIter]=spin1;
		}
	      tPartIter++;	//next particle
	      for(tIter2=0;tIter2<50;tIter2++) str[tIter2]='\0';
	    }
	}
      //END OF DATA
    }
  mParticleCount=tPartIter;
  //////  END OF "TABLES.M" /////////


  //////	START OF "i200STAR.m"//////

  double izospin1,izospin2;	//help to calculate izospin, if niecalkowity
  //    input=fopen("i200STAR.m","r+");
  ifstream in1("i200STAR.m");

  for(tIter2=0;tIter2<369;tIter2++)
    for(int jj=0;jj<2;jj++) mDecayBuffer[tIter2][jj]=0; 
  
  for(;;)
    {
      j=3; //first letter of particle in line in file
      for(tIter=0;tIter<50;tIter++) str[tIter]='\0';
      for(tIter=0;tIter<20;tIter++) str1[tIter]='\0';
      tIter=0;
	
      in1.getline(str,200);
      
      if(str[0] == 'x') break;

      if(str[0] == 'f' && str[1] == 'i')
	{	
	  //name
	  for(;;)	
	    {
	      if(str[j] == ',') 
		{
		  j++;	//to skip ","
		  break;
		}
	      str1[tIter++]=str[j++];
	    }
	  for(tIter=0;tIter<369;tIter++)
	    {	    
	      if(check(str1,mNameBuffer[tIter]))
		{
		  //barion number
		  for(tIter2=0;tIter2<20;tIter2++) str1[tIter2]='\0';
		  if(str[j] == '-')
	    	    { 
		      *str1=str[j+1];
		      mBarionBuffer[tIter]=atoi(str1);
		      mBarionBuffer[tIter]=-mBarionBuffer[tIter];
		    }
		  if(str[j] != '-') 
	    	    {
		      *str1=str[j];
		      mBarionBuffer[tIter]=atoi(str1);
	    	    }
		  if(str[j] == '-') j+=3;
		  else j+=2;

		  //strange number
		  for(tIter2=0;tIter2<20;tIter2++) str1[tIter2]='\0';
		  if(str[j] == '-')
	    	    {
		      *str1=str[j+1];
		      mStrangeBuffer[tIter]=atoi(str1);
		      mStrangeBuffer[tIter]=-mStrangeBuffer[tIter];
		    }
		  if(str[j] != '-') 
		    {
		      *str1=str[j];
		      mStrangeBuffer[tIter]=atoi(str1);
		    }
		  if(str[j] == '-') j+=3;
		  else j+=2;

		  //izospin3 number
		  for(tIter2=0;tIter2<20;tIter2++) str1[tIter2]='\0';
		  if(str[j+1] == '/' && str[j+3] == ']') 
	    	    {
		      *str1=str[j];
		      izospin1=atoi(str1);
		      *str1=str[j+2];
		      izospin2=atoi(str1);
		      mI3Buffer[tIter]=izospin1/izospin2;
		    }
		  if(str[j] == '-' && str[j+2] == '/' && str[j+4] == ']') 
		    {
		      *str1=str[j+1];
		      izospin1=atoi(str1);
		      *str1=str[j+3];
		      izospin2=atoi(str1);
		      mI3Buffer[tIter]=-izospin1/izospin2;
		    }
		  if(str[j] == '-' && str[j+2] == ']')
		    {
		      *str1=str[j+1];
		      izospin1=atoi(str1);
		      mI3Buffer[tIter]=-izospin1;
		    }
		  if(str[j+1] == ']') 
		    {
		      *str1=str[j];
		      mI3Buffer[tIter]=atof(str1);
		    }
		  break;
		}
	    }
	}
	
      //DECAY CHANNELS
	
      tIter=0;
      //TWO-BODY DECAYS
      if(str[0] == 's' && str[1] == 'e' && str[2] == '[')
	{
	  // Reading in the decay channels
	  char *tLBrackert, *tFirstComma, *tSecondComma, *tThirdComma, *tRBracket;
	  
	  tLBrackert = strchr(str,'[');
	  tFirstComma = strchr(str,',');
	  tSecondComma = strchr(tFirstComma+1,',');
	  tThirdComma = strchr(tSecondComma+1,',');
	  tRBracket = strchr(tThirdComma,']');

	  if (!(tLBrackert && tFirstComma && tSecondComma && tThirdComma && tRBracket))
	    PRINT_DEBUG_1("Malformed line!: " << str);
	  
	  char *tFather = new char[tFirstComma-tLBrackert];
	  strncpy(tFather, tLBrackert+1,   tFirstComma-tLBrackert-1);
	  char *tDaughter1 = new char[tSecondComma-tFirstComma];
	  strncpy(tDaughter1, tFirstComma+1,  tSecondComma-tFirstComma-1);
	  char *tDaughter2 = new char[tThirdComma-tSecondComma];
	  strncpy(tDaughter2, tSecondComma+1, tThirdComma-tSecondComma-1);
	  char *tBRatio = new char[tRBracket-tThirdComma];
	  strncpy(tBRatio, tThirdComma+1,  tRBracket-tThirdComma-1);
	  
	  // Getting the ratio
	  char *tMiddle, *tRatComponent;
	  double tRatio = 1.0;
	  
	  tMiddle = strchr(tBRatio,'_');
	  if (!tMiddle) tMiddle = strchr(tBRatio,' ');
	  while(tMiddle)
	    {
	      tRatComponent = new char[tMiddle-tBRatio+1];
	      strncpy(tRatComponent, tBRatio, tMiddle-tBRatio);
	      if (strchr(tRatComponent,'/'))
		tRatio *= atof(tRatComponent)/atof(tRatComponent+2);
	      else
		tRatio *= atof(tRatComponent);
	      tBRatio = tMiddle+1;		
	      tMiddle = strchr(tBRatio,'_');
	      if (!tMiddle) tMiddle = strchr(tBRatio,' ');
	      delete [] tRatComponent;
	    }
	  if (strchr(tBRatio,'/'))
	    tRatio *= atof(tBRatio)/atof(tBRatio+2);
	  else
	    tRatio *= atof(tBRatio);

	  delete [] tFather;
	  delete [] tDaughter1;
	  delete [] tDaughter2;
	  delete [] tBRatio;
	}

      //THREE-BODY DECAYS
      j++;	//because se3[]
      if(str[0] == 's' && str[1] == 'e' && str[2] == '3' && str[3] == '[')
	{
	  for(;;)	
	    {
	      if(str[j] == ',') 
		{
		  j++;	//to skip ","
		  break;
		}
	      str1[tIter++]=str[j++];	//name
	    }

	  for(tIter=0;tIter<369;tIter++)
	    {	    
	      if(check(str1,mNameBuffer[tIter]))
		{
		  mDecayBuffer[tIter][1]++;
		  break;
		}
	    }
	}


    }
  //////	END OF "i200STAR.m"//////
    
  ParticleType *tPartBuf;
  int tNum, tPart;
  
  for(tPart=0;tPart<mParticleCount;tPart++)
    {
      tPartBuf = new ParticleType();
      tPartBuf->SetName(mNameBuffer[tPart]);
      tPartBuf->SetMass(mMassBuffer[tPart]);
      tPartBuf->SetGamma(mGammaBuffer[tPart]);
      tPartBuf->SetSpin(mSpinBuffer[tPart]);
      tPartBuf->SetBarionN(mBarionBuffer[tPart]);
      tPartBuf->SetStrangeness(mStrangeBuffer[tPart]);
      tPartBuf->SetI3(mI3Buffer[tPart]);
      tPartBuf->SetDecayChannelCount2(mDecayBuffer[tPart][0]);
      tPartBuf->SetDecayChannelCount3(mDecayBuffer[tPart][1]);
      tPartBuf->SetNumber(tPart);
      tNum = mDB->AddParticleType(tPartBuf);
    }

  ifstream in2("i200STAR.m");
  while (in2.getline(str,200))
    {
      tIter=0;
      //TWO-BODY DECAYS
      if(str[0] == 's' && str[1] == 'e' && str[2] == '[')
	{
	  // Reading in the decay channels
	  char *tLBrackert, *tFirstComma, *tSecondComma, *tThirdComma, *tRBracket;
	  
	  tLBrackert = strchr(str,'[');
	  tFirstComma = strchr(str,',');
	  tSecondComma = strchr(tFirstComma+1,',');
	  tThirdComma = strchr(tSecondComma+1,',');
	  tRBracket = strchr(tThirdComma,']');
	  if (!(tLBrackert && tFirstComma && tSecondComma && tThirdComma && tRBracket))
	    PRINT_DEBUG_1("Malformed line!: " << str);
	  
	  char *tFather = new char[tFirstComma-tLBrackert];
	  strncpy(tFather, tLBrackert+1,   tFirstComma-tLBrackert-1);
	  char *tDaughter1 = new char[tSecondComma-tFirstComma];
	  strncpy(tDaughter1, tFirstComma+1,  tSecondComma-tFirstComma-1);
	  char *tDaughter2 = new char[tThirdComma-tSecondComma];
	  strncpy(tDaughter2, tSecondComma+1, tThirdComma-tSecondComma-1);
	  char *tBRatio = new char[tRBracket-tThirdComma];
	  strncpy(tBRatio, tThirdComma+1,  tRBracket-tThirdComma-1);
	  
	  // Getting the ratio
	  char *tMiddle, *tRatComponent;
	  double tRatio = 1.0;
	  
	  tMiddle = strchr(tBRatio,'_');
	  if (!tMiddle) tMiddle = strchr(tBRatio,' ');
	  while(tMiddle)
	    {
	      tRatComponent = new char[tMiddle-tBRatio];
	      strncpy(tRatComponent, tBRatio, tMiddle-tBRatio);
	      if (strchr(tRatComponent,'/'))
		tRatio *= atof(tRatComponent)/atof(tRatComponent+2);
	      else
		tRatio *= atof(tRatComponent);
	      tBRatio = tMiddle+1;		
	      tMiddle = strchr(tBRatio,'_');
	      if (!tMiddle) tMiddle = strchr(tBRatio,' ');
	      delete [] tRatComponent;
	    }
	  if (strchr(tBRatio,'/'))
	    tRatio *= atof(tBRatio)/atof(tBRatio+2);
	  else
	    tRatio *= atof(tBRatio);
	  
	  DecayChannel *newChannel = new DecayChannel(tRatio, mDB->GetParticleTypeIndex(tDaughter1), mDB->GetParticleTypeIndex(tDaughter2), -1);
	  // 	  if (mDB->GetParticleType(tDaughter1)->GetMass() + 
	  // 	      mDB->GetParticleType(tDaughter2)->GetMass()
	  // 	      < mDB->GetParticleType(tFather)->GetMass()) {
	  // 	    (mDB->GetParticleType(tFather))->AddDecayChannel(*newChannel);
	  // 	  }
	  (mDB->GetParticleType(tFather))->AddDecayChannel(*newChannel);

	  delete newChannel;
	}
      
      if(str[0] == 's' && str[1] == 'e' && str[2] == '3')
	{
	  // Reading in the decay channels
	  char *tLBrackert, *tFirstComma, *tSecondComma, *tThirdComma, *tFourthComma, *tRBracket;
	  
	  tLBrackert = strchr(str,'[');
	  tFirstComma = strchr(str,',');
	  tSecondComma = strchr(tFirstComma+1,',');
	  tThirdComma = strchr(tSecondComma+1,',');
	  tFourthComma = strchr(tThirdComma+1,',');
	  tRBracket = strchr(tThirdComma,']');

	  if (!(tLBrackert && tFirstComma && tSecondComma && tThirdComma && tFourthComma && tRBracket))
	    PRINT_DEBUG_1("Malformed line!: " << str);
	  
	  char *tFather = new char[tFirstComma-tLBrackert];
	  strncpy(tFather, tLBrackert+1,   tFirstComma-tLBrackert-1);
	  char *tDaughter1 = new char[tSecondComma-tFirstComma];
	  strncpy(tDaughter1, tFirstComma+1,  tSecondComma-tFirstComma-1);
	  char *tDaughter2 = new char[tThirdComma-tSecondComma];
	  strncpy(tDaughter2, tSecondComma+1, tThirdComma-tSecondComma-1);
	  char *tDaughter3 = new char[tFourthComma-tThirdComma];
	  strncpy(tDaughter3, tThirdComma+1,  tFourthComma-tThirdComma-1);
	  char *tBRatio = new char[tRBracket-tFourthComma];
	  strncpy(tBRatio, tFourthComma+1, tRBracket-tFourthComma-1);
	  
	  // Getting the ratio
	  char *tMiddle, *tRatComponent;
	  double tRatio = 1.0;
	  
	  tMiddle = strchr(tBRatio,'_');
	  if (!tMiddle) tMiddle = strchr(tBRatio,' ');
	  while(tMiddle)
	    {
	      tRatComponent = new char[tMiddle-tBRatio+1];
	      strncpy(tRatComponent, tBRatio, tMiddle-tBRatio);
	      if (strchr(tRatComponent,'/'))
		tRatio *= atof(tRatComponent)/atof(tRatComponent+2);
	      else
		tRatio *= atof(tRatComponent);
	      tBRatio = tMiddle+1;		
	      tMiddle = strchr(tBRatio,'_');
	      if (!tMiddle) tMiddle = strchr(tBRatio,' ');
	      delete [] tRatComponent;
	    }
	  if (strchr(tBRatio,'/'))
	    tRatio *= atof(tBRatio)/atof(tBRatio+2);
	  else
	    tRatio *= atof(tBRatio);
	  
	  DecayChannel *newChannel = new DecayChannel(tRatio, mDB->GetParticleTypeIndex(tDaughter1), mDB->GetParticleTypeIndex(tDaughter2), mDB->GetParticleTypeIndex(tDaughter3));
	  // 	  if (mDB->GetParticleType(tDaughter1)->GetMass() + 
	  // 	      mDB->GetParticleType(tDaughter2)->GetMass() +
	  // 	      mDB->GetParticleType(tDaughter3)->GetMass()
	  // 	      < mDB->GetParticleType(tFather)->GetMass())
	  (mDB->GetParticleType(tFather))->AddDecayChannel(*newChannel);

	  delete tFather;
	  delete tDaughter1;
	  delete tDaughter2;
	  delete tDaughter3;
	  delete tBRatio;
	}
    }

  ifstream in3("pdgcodes.m");
  while (in3.getline(str,200))
    {
      string tName;
      int tCode;
      
      std::stringstream tIS(str);
      tIS >> tName >> tCode;
      mDB->GetParticleType(tName)->SetPDGCode(tCode);
    }
  
}


void Parser::ReadShare()
{
  char str[50];
  char str1[200];
    
  ParticleType *tPartBuf;
  int tNum;
           
  PRINT_DEBUG_1("Reading from |"<<(mInputDir+"/"+"particles.data").Data()<<"|");
  ifstream in((mInputDir+"/"+"particles.data").Data());
    
  int number=0;
  if ((in) && (in.is_open()))
    {
      //START OF HEAD-LINE
      in.ignore(200,'\n');
      in.ignore(200,'\n');
      in.ignore(200,'\n');
      //END OF HEAD-LINE
      
      while (in>>str)
	{
	  if (/*(*str == '#')||*/(*str<65)||(*str>122))
	    {
	      in.getline(str1,200);
	      continue;
	    }
	  double mass, gamma, spin, I3, I, q, s, aq, as, c, ac, mc;
	  
	  in>>mass>>gamma>>spin>>I>>I3>>q>>s>>aq>>as>>c>>ac>>mc;
	  number++;
	  PRINT_DEBUG_2(number<<" "<<str<<" "<<mass<<" "<<gamma<<" "<<spin<<" "<<I<<" "<<I3<<" "<<q<<" "<<aq<<" "<<s<<" "<<as<<" "<<c<<" "<<ac<<" "<<mc);
	  
	  tPartBuf = new ParticleType();
	  tPartBuf->SetName(str);
	  tPartBuf->SetMass(mass);
	  tPartBuf->SetGamma(gamma);
	  tPartBuf->SetSpin(spin);
	  tPartBuf->SetBarionN((int) ((q+s+c)/3. - (aq+as+ac)/3.) );
	  tPartBuf->SetCharmN((int) (c - ac));
	  tPartBuf->SetStrangeness((int) (as-s));
	  tPartBuf->SetI(I);
	  tPartBuf->SetI3(I3);
	  tPartBuf->SetPDGCode((int) mc);
	  tPartBuf->SetNumber(number);
	  tNum = mDB->AddParticleType(tPartBuf);
	}
      in.close();
    }
  else 
    {
      PRINT_MESSAGE("File "<<(mInputDir+"/"+"particles.data").Data()<<" containing particle data not found!");
      PRINT_MESSAGE("Please set the correct path to this file in the input parameter file");
      PRINT_MESSAGE("Aborting!");
      exit(0);
    }
  
  ifstream in2((mInputDir+"/decays.data").Data());
  if ((in2) && (in2.is_open()))
    {
      //START OF HEAD-LINE
      in2.ignore(200,'\n');
      in2.ignore(200,'\n');
      in2.ignore(200,'\n');
      //END OF HEAD-LINE

      char tFather[50], tDaughter1[50], tDaughter2[50], tDaughter3[50];
      double tBRatio, tRatio;
      int CGcoeff; // complete branching ratio by Clebsch-Gordan coefficient: 0-no 1-yes
	
      while (in2>>str)
	{
	  if (*str == '#')
	    {
	      in2.getline(str1,200);
	      continue;
	    }
	  in2>>tDaughter1>>tDaughter2>>tDaughter3;
	  if (!mDB->ExistsParticleType(tDaughter1)) {
	    PRINT_MESSAGE("Did not find the daughter 1 particle: " << tDaughter1);
	    PRINT_MESSAGE("Not adding channel");
	    in2.getline(str1,200);
	    continue;
	  }
	  if (!mDB->ExistsParticleType(tDaughter2)) {
	    PRINT_MESSAGE("Did not find the daughter 2 particle: " << tDaughter2);
	    PRINT_MESSAGE("Not adding channel");
	    in2.getline(str1,200);
	    continue;
	  }
	  if ((*tDaughter3>65)&&(*tDaughter3<122)&&(!mDB->ExistsParticleType(tDaughter3))) {
	    PRINT_MESSAGE("Did not find the daughter 3 particle: " << tDaughter3);
	    PRINT_MESSAGE("Not adding channel");
	    in2.getline(str1,200);
	    continue;
	  }

	  strcpy(tFather,str);
	  PRINT_DEBUG_2(tFather<<"\t"<<tDaughter1<<"\t"<<tDaughter2<<"\t");
	  if ((*tDaughter3>65)&&(*tDaughter3<122)) // check if first char is a letter - if yes then 3-body decay
	    {
	      in2>>tBRatio>>CGcoeff;
	      PRINT_DEBUG_2(tDaughter3<<" (3-body decay)\t");
	      if (mDB->ExistsParticleType(tFather)) {
		mDB->GetParticleType(tFather)->SetDecayChannelCount3(mDB->GetParticleType(tFather)->GetDecayChannelCount3()+1);
	      
		tRatio=tBRatio;
		DecayChannel *newChannel = new DecayChannel(tRatio, mDB->GetParticleTypeIndex(tDaughter1), mDB->GetParticleTypeIndex(tDaughter2), mDB->GetParticleTypeIndex(tDaughter3));
		if (mDB->GetParticleType(tDaughter1)->GetMass() + 
		    mDB->GetParticleType(tDaughter2)->GetMass() +
		    mDB->GetParticleType(tDaughter3)->GetMass()
		    < mDB->GetParticleType(tFather)->GetMass())
		  (mDB->GetParticleType(tFather))->AddDecayChannel(*newChannel);

		delete newChannel;
		
		tRatio=tBRatio;
		PRINT_DEBUG_2(tBRatio << '\t' << tRatio);
	      }
	      else {
		PRINT_MESSAGE("Did not find the father particle: " << tFather);
		PRINT_MESSAGE("Not adding channel");
	      }
	    }
	  else // 2-body decay	    
	    {
	      tBRatio=atof(tDaughter3);
	    
	      in2>>CGcoeff;
	      PRINT_DEBUG_2(" (2-body decay)\t");
	      if (mDB->ExistsParticleType(tFather)) {
		mDB->GetParticleType(tFather)->SetDecayChannelCount2(mDB->GetParticleType(tFather)->GetDecayChannelCount2()+1);
		
		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if (CGcoeff) // complete branching ratio by Clebsch-Gordan coefficient
		  {
		    double j1, m1, j2, m2, J, M, CB;
		    J=mDB->GetParticleType(tFather)->GetI();
		    M=mDB->GetParticleType(tFather)->GetI3();
		    j1=mDB->GetParticleType(tDaughter1)->GetI();
		    m1=mDB->GetParticleType(tDaughter1)->GetI3();
		    j2=mDB->GetParticleType(tDaughter2)->GetI();
		    m2=mDB->GetParticleType(tDaughter2)->GetI3();
		    PRINT_DEBUG_2(" "<<J<<" "<<M<<" "<<j1<<" "<<m1<<" "<<j2<<" "<<m2);
		    
		    CB = ClebschGordan(J, M, j1, m1, j2, m2);
		    tRatio = CB*CB * tBRatio;
		    
		    // Multiply the Clebsh by two?
		    // The same spin, mass, strangeness
		    // and different I3?
		    if ((fabs(mDB->GetParticleType(tDaughter1)->GetSpin() - mDB->GetParticleType(tDaughter2)->GetSpin()) < 0.01) &&
			(fabs(mDB->GetParticleType(tDaughter1)->GetMass() - mDB->GetParticleType(tDaughter2)->GetMass()) < 0.01) &&
			(mDB->GetParticleType(tDaughter1)->GetStrangeness() == mDB->GetParticleType(tDaughter2)->GetStrangeness()) &&
			(fabs(mDB->GetParticleType(tDaughter1)->GetI3() - mDB->GetParticleType(tDaughter2)->GetI3()) > 0.01))
		      {
			PRINT_DEBUG_2("Multuplying Clebsch by two for " << tFather << "->" << tDaughter1 << "+" << tDaughter2);
			tRatio *= 2.0;
		      }
		    
		    PRINT_DEBUG_2(CB << '\t' << tBRatio << '\t' << tRatio<<"\t"<<CGcoeff);
		  }
		
		else
		  {
		    tRatio=tBRatio;
		    PRINT_DEBUG_2(tBRatio << '\t' << tRatio);
		  }
		DecayChannel *newChannel = new DecayChannel(tRatio, mDB->GetParticleTypeIndex(tDaughter1), mDB->GetParticleTypeIndex(tDaughter2), -1);
		if (mDB->GetParticleType(tDaughter1)->GetMass() + 
		    mDB->GetParticleType(tDaughter2)->GetMass()
		    < mDB->GetParticleType(tFather)->GetMass()) 
		  {
		    (mDB->GetParticleType(tFather))->AddDecayChannel(*newChannel);
		    PRINT_DEBUG_2("Added channel " << newChannel << " " << mDB->GetParticleTypeIndex(tFather) << " " << mDB->GetParticleTypeIndex(tDaughter1) << " " << mDB->GetParticleTypeIndex(tDaughter2));
		  }
		else 
		  {
		    
		    PRINT_DEBUG_2("Masses do not match! Not adding channel " << newChannel);

		    delete newChannel;
		  }
	      }
	      else {
		PRINT_MESSAGE("Did not find the father particle: " << tFather);
		PRINT_MESSAGE("Not adding channel");
	      }
	    }
	}
      in2.close();
    }
  else {
    PRINT_MESSAGE("File "<<(mInputDir+"/decays.data").Data()<<" with particle decay channels not found!");
    PRINT_MESSAGE("No particle decays will be simulated");
  }
}

double 
Parser::ClebschGordan(double aJot, double aEm, double aJot1, double aEm1, double aJot2, double aEm2)
{
  int mint, maxt;
  double cgc = 0.0;
  int titer;
  double coef;

  maxt = lrint(aJot1 + aJot2 - aJot);
  mint = 0;
  if (lrint(aJot1 - aEm1) < maxt) maxt = lrint(aJot1 - aEm1);
  if (lrint(aJot2 + aEm2) < maxt) maxt = lrint(aJot2 + aEm2);
  if (lrint(-(aJot-aJot2+aEm1)) > mint) mint = lrint(-(aJot-aJot2+aEm1));
  if (lrint(-(aJot-aJot1-aEm2)) > mint) mint = lrint(-(aJot-aJot1-aEm2));

  PRINT_DEBUG_3("mint " << mint << " " <<  aJot1 << " " << aEm1);
  PRINT_DEBUG_3("maxt " << maxt << " " <<  aJot2 << " " << aEm2);

  for (titer = mint; titer<=maxt; titer ++)
    {
      coef = TMath::Power(-1, titer);
      PRINT_DEBUG_3("coef1 " << coef); 
      coef *= TMath::Sqrt((2*aJot+1)*
			  factorials[lrint(aJot1+aEm1)] *
			  factorials[lrint(aJot1-aEm1)] *
			  factorials[lrint(aJot2+aEm2)] *
			  factorials[lrint(aJot2-aEm2)] *
			  factorials[lrint(aJot+aEm)] *
			  factorials[lrint(aJot-aEm)]);
      PRINT_DEBUG_3("coef2 " << coef); 
      coef /= (factorials[titer] *
	       factorials[lrint(aJot1+aJot2-aJot-titer)] *
	       factorials[lrint(aJot1-aEm1-titer)] *
	       factorials[lrint(aJot2+aEm2-titer)] *
	       factorials[lrint(aJot-aJot2+aEm1+titer)] *
	       factorials[lrint(aJot-aJot1-aEm2+titer)]);
      PRINT_DEBUG_3("coef3 " << coef); 
      
      cgc += coef;
    }

  cgc *= DeltaJ(aJot1, aJot2, aJot);

  return cgc;
}

double 
Parser::DeltaJ(double aJot1, double aJot2, double aJot)
{
  double res = TMath::Sqrt(1.0 * 
			   factorials[lrint(aJot1+aJot2-aJot)] * 
			   factorials[lrint(aJot1-aJot2+aJot)] * 
			   factorials[lrint(-aJot1+aJot2+aJot)] / 
			   factorials[lrint(aJot1+aJot2+aJot+1)]);
  
  return res;
}

void   
Parser::ReadParameters()
{
  // Read the input directory
  try {
    mInputDir = sRPInstance->getPar("InputDirSHARE");
  }
  catch (STR tError) {
    PRINT_DEBUG_1("Parser::ReadParameters - Caught exception " << tError);
    PRINT_MESSAGE("Did not find SHARE input file location.");
    PRINT_MESSAGE("Using default: '../share'");
    mInputDir = "../share";
  }
}
