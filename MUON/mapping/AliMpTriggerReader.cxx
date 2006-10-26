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
* about the suitability of this software for any purpeateose. It is      *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id$
// $MpId: AliMpTriggerReader.cxx,v 1.4 2006/05/24 13:58:52 ivana Exp $

#include "AliMpTriggerReader.h"
#include "AliMpMotifReader.h"
#include "AliMpFiles.h"
#include "AliMpMotifType.h"
#include "AliMpPCB.h"
#include "AliMpSlat.h"
#include "AliMpMotifMap.h"
#include "AliMpMotifPosition.h"
#include "AliMpMotif.h"
#include "AliMpHelper.h"
#include "AliMpSt345Reader.h"
#include "AliMpTrigger.h"

#include "AliLog.h"

#include "Riostream.h"
#include "TClass.h"
#include "TObjString.h"
#include "TList.h"
#include "TString.h"

#include <sstream>
#include <cassert>

/// 
/// \class AliMpTriggerReader
/// Read trigger slat ASCII files
/// Basically provides 2 static methods:
/// - AliMpTrigger* ReadSlat()
/// - AliMpPCB* ReadPCB()
///
/// \author Laurent Aphecetche

/// \cond CLASSIMP
ClassImp(AliMpTriggerReader)
/// \endcond

TMap AliMpTriggerReader::fgPCBMap;
TMap AliMpTriggerReader::fgLocalBoardMap;

const TString AliMpTriggerReader::fgkKeywordLayer("LAYER");
const TString AliMpTriggerReader::fgkKeywordScale("SCALE");
const TString AliMpTriggerReader::fgkKeywordPcb("PCB");  
const TString AliMpTriggerReader::fgkKeywordFlipX("FLIP_X");
const TString AliMpTriggerReader::fgkKeywordFlipY("FLIP_Y");

//_____________________________________________________________________________
AliMpTriggerReader::AliMpTriggerReader() : TObject()
{
  //
  // Default ctor.
  //
} 

//_____________________________________________________________________________
AliMpTriggerReader::~AliMpTriggerReader()
{
  //
  // Dtor.
  //
  fgPCBMap.Delete();
  fgLocalBoardMap.Delete();
}

//_____________________________________________________________________________
AliMpSlat*
AliMpTriggerReader::BuildSlat(const char* slatName,
                              AliMpPlaneType planeType,
                              const TList& lines,
                              Double_t scale)
{
  // Construct a slat from the list of lines, taking into account
  // the scale factor

  AliMpSlat* slat = new AliMpSlat(slatName, planeType);
    
  TIter it(&lines);
  TObjString* osline;
  while ( ( osline = (TObjString*)it.Next() ) )
  {
    // note that at this stage lines should not be empty.
    TString sline(osline->String());
    
    TObjArray* tokens = sline.Tokenize(' ');
    
    TString& keyword = ((TObjString*)tokens->At(0))->String();
    
    if ( keyword == fgkKeywordPcb )
    {
      if ( tokens->GetEntriesFast() != 3 )
      {
        AliErrorClass(Form("Syntax error : expecting PCB type localboard-list"
                           " in following line:\n%s",sline.Data()));
        delete slat;
        delete tokens;
        return 0;
      }
      TString pcbName = ((TObjString*)tokens->At(1))->String();
      
      TObjArray* localBoardList = ((TObjString*)tokens->At(2))->String().Tokenize(',');
      
      if ( scale != 1.0 )
      {
        std::ostringstream s;
        s << pcbName.Data() << "x" << scale;
        pcbName = s.str().c_str();
      }
      
      AliMpPCB* pcbType = PCB(pcbName.Data());	  
      if (!pcbType)
      {
        AliErrorClass(Form("Cannot read pcbType=%s",pcbName.Data()));
        delete slat;
        return 0;
      }      

      TArrayI allLocalBoards;
      
      for ( Int_t ilb = 0; ilb < localBoardList->GetEntriesFast(); ++ilb)
      {
        TArrayI localBoardNumbers;
        TString& localBoards = ((TObjString*)localBoardList->At(ilb))->String();
        Ssiz_t pos = localBoards.First('-');
        if ( pos < 0 ) 
        {
          pos = localBoards.Length();
        }
        AliMpHelper::DecodeName(localBoards(pos-1,localBoards.Length()-pos+1).Data(),
                                ';',localBoardNumbers);      
        for ( int i = 0; i < localBoardNumbers.GetSize(); ++i )
        {
          std::ostringstream name;
          name << localBoards(0,pos-1) << localBoardNumbers[i];
          AliDebugClass(3,name.str().c_str());
          localBoardNumbers[i] = LocalBoardNumber(name.str().c_str());
          AliDebugClass(3,Form("LOCALBOARDNUMBER %d\n",localBoardNumbers[i]));
          allLocalBoards.Set(allLocalBoards.GetSize()+1);
          allLocalBoards[allLocalBoards.GetSize()-1] = localBoardNumbers[i];
          if (localBoardNumbers[i] < 0 )
          {
            AliErrorClass(Form("Got a negative local board number in %s ? Unlikely"
                               " to be correct... : %s\n",slatName,name.str().c_str()));
          }
        }
      }
      delete tokens;
      delete localBoardList;
      slat->Add(pcbType,allLocalBoards);
    }
  }
  
  if ( slat->DX()== 0 || slat->DY() == 0 )
  {
    AliFatalClass(Form("Slat %s has invalid null size\n",slat->GetID()));
  }
  return slat;
}

//_____________________________________________________________________________
TString
AliMpTriggerReader::GetBoardNameFromPCBLine(const TString& s)
{
  // Decode the string to get the board name
  TString boardName;
  
  TObjArray* tokens = s.Tokenize(' ');
  
  TString& keyword = ((TObjString*)tokens->At(0))->String();

  if ( keyword == fgkKeywordPcb &&
       tokens->GetEntriesFast() == 3 )
  {
    boardName = ((TObjString*)tokens->At(2))->String();
  }
  
  delete tokens;
  
  return boardName;
}
  
//_____________________________________________________________________________
void
AliMpTriggerReader::FlipLines(TList& lines, Bool_t flipX, Bool_t flipY,
                              Int_t srcLine, Int_t destLine)
{
  //
  // Change the local board names contained in lines, 
  // to go from right to left, and/or
  // from top to bottom
  //
 

  if ( flipX )
  {
    // Simply swaps R(ight) and L(eft) in the first character of 
    // local board names

    TObjString* oline;
    TIter it(&lines);
    while ( ( oline = (TObjString*)it.Next() ) )
    {
      TString& s = oline->String();
      if ( s.Contains("RC") ) 
      {
        // Change right to left
        s.ReplaceAll("RC","LC");
      }
      else if ( s.Contains("LC") )
      {
        // Change left to right
        s.ReplaceAll("LC","RC");
      }
    }
  }
  
  if ( flipY )
  {
    // Change line number, according to parameters srcLine and destLine
    // Note that because of road opening (for planes 3 and 4 at least),
    // we loop for srcLine +-1
    //
    for ( Int_t line = -1; line <=1; ++line )
    {
      std::ostringstream src,dest;
      src << "L" << srcLine+line;
      dest << "L" << destLine-line;
      if ( src.str() == dest.str() ) continue;
      
      for ( Int_t i = 0; i < lines.GetSize(); ++i )
      {
        TObjString* oline = (TObjString*)lines.At(i);
        
        TString& s = oline->String();
        
        if ( !s.Contains(fgkKeywordPcb) )
        {
          // Only consider PCB lines.
          continue;
        }
        
        if ( s.Contains(src.str().c_str()) )
        {
          AliDebugClass(4,Form("Replacing %s by %s in %s\n",
                               src.str().c_str(),dest.str().c_str(),s.Data()));
          
          s.ReplaceAll(src.str().c_str(),dest.str().c_str());
          
          AliDebugClass(4,s.Data());
          
          TString boardName(GetBoardNameFromPCBLine(s));
          
          if ( line )
          {
            // We must also change board numbers, with the tricky
            // thing that up and down must be swapped...
            // Up can only be 1 card so it must be B1
            // Down must be the uppper card of the line before, so
            // the biggest possible board number for this Line,Column
            
            if (line>0)
            {
                // force to B1
              AliDebugClass(4,Form("Forcing B1 in %s\n",s.Data()));
              s.ReplaceAll(boardName(boardName.Length()-2,2),"B1");
              AliDebugClass(4,s.Data());
            }
            else
            {
              // find the largest valid board number
              for ( int b = 4; b>=1; --b )
              {
                std::ostringstream bs;
                bs << boardName(0,boardName.Length()-1) << b;
                if ( LocalBoardNumber(bs.str().c_str()) >= 0 )
                {
                  AliDebugClass(4,Form("Replacing %s by %s in %s\n",
                                  boardName(boardName.Length()-2,2).Data(),
                                  Form("B%d",b),
                                  s.Data()));
                  s.ReplaceAll(boardName(boardName.Length()-2,2),
                               Form("B%d",b));
                  AliDebugClass(4,s);
                  break;
                }
              }
            }  
            // Check that the replacement we did is ok. If not,
            // skip the line.
            Int_t lbn = LocalBoardNumber(GetBoardNameFromPCBLine(s));
            if ( lbn < 0 )
            {
              AliDebugClass(4,Form("Removing line %s\n",s.Data()));
              lines.Remove(oline);
            }
            
          } // if (line)          
        }
      }    
    }
  }
}

//___________________________________________________________________________
Int_t
AliMpTriggerReader::IsLayerLine(const TString& sline)
{
  // Whether sline contains LAYER keyword

  if ( sline.BeginsWith(fgkKeywordLayer) )
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

//___________________________________________________________________________
Int_t
AliMpTriggerReader::DecodeFlipLine(const TString& sline,
                                   TString& slatType2,
                                   Bool_t& flipX, Bool_t& flipY)
{
  // Decode a line containing FLIP_X and/or FLIP_Y keywords

  Ssiz_t blankPos = sline.First(' ');
  if ( blankPos < 0 ) return 0;
  
  TString keyword(sline(0,blankPos));
  
  if ( keyword == fgkKeywordFlipX )
  {
    flipX = kTRUE;
  } else if ( keyword == fgkKeywordFlipY )
  {
    flipY = kTRUE;
  }
  else
  {
    return 0;
  }
  
  slatType2 = sline(blankPos+1,sline.Length()-blankPos-1);
  return 1;
}

//___________________________________________________________________________
Int_t
AliMpTriggerReader::DecodeScaleLine(const TString& sline, 
                                    Double_t& scale, TString& slatType)
{
  // Decode sline containing SCALE keyword

  if ( sline(0,fgkKeywordScale.Length()) == fgkKeywordScale )
  {
    TString tmp(sline(fgkKeywordScale.Length()+1,
                      sline.Length()-fgkKeywordScale.Length()-1));
    Ssiz_t blankPos = tmp.First(' ');
    if ( blankPos < 0 )
    {
      AliErrorClass(Form("Syntax error in slat file, should get a slatType after "
                    " SCALE keyword : %s\n",tmp.Data()));
      return -1;
    }
    else
    {
      slatType = tmp(0,blankPos);
      scale = TString(tmp(blankPos+1,tmp.Length()-blankPos-1)).Atof();
      return 1;
    }
  }
  scale = 1.0;
  return 0;
}

//_____________________________________________________________________________
Int_t
AliMpTriggerReader::GetLine(const TString& slatType)
{
  //
  // Assuming slatType is a 4 character string of the form XSLN
  // where X=1,2,3 or 4
  // S = R or L
  // N is the line number
  // returns N
  if ( isdigit(slatType[0]) && 
       ( slatType[1] == 'R' || slatType[1] == 'L' ) &&
       slatType[2] == 'L' )
  {
    return atoi(slatType(3,1).Data());
  }
  return -1;
}

//_____________________________________________________________________________
int
AliMpTriggerReader::LocalBoardNumber(const char* localBoardName)
{
  // From local board name to local board number

  if ( !fgLocalBoardMap.GetSize() ) 
  {
    ReadLocalBoardMapping();
  }
  
  TPair* pair = (TPair*)fgLocalBoardMap.FindObject(localBoardName);
  
  if (pair)
  {
    return atoi(((TObjString*)pair->Value())->String().Data());
  }
  return -1;
}

//_____________________________________________________________________________
AliMpPCB*
AliMpTriggerReader::PCB(const char* pcbType)
{
  //
  // Get access to an AliMpPCB object, given its type (e.g. N1, SB2, etc...)
  //
  // Note that the returned object is either a new one (read from file) or a 
  // reused one if it is already present in the internal map.
  //
  
  TPair* pair = (TPair*)fgPCBMap.FindObject(pcbType);
  if ( pair )
  {
    AliDebugClass(2,Form("Getting pcb %s from internal map",pcbType));
    return (AliMpPCB*)pair->Value();
  }
  else
  {
    AliDebugClass(2,Form("Reading pcb %s from file",pcbType));
    return ReadPCB(pcbType);
  }
}

//_____________________________________________________________________________
void 
AliMpTriggerReader::ReadLines(const char* slatType,
                              AliMpPlaneType planeType,
                              TList& lines,
                              Double_t& scale,
                              Bool_t& flipX, Bool_t& flipY,
                              Int_t& srcLine, Int_t& destLine)
{
  //
  // Reads in lines from file for a given slat
  // Returns the list of lines (lines), together with some global
  // information as the scale, whether to flip the lines, etc...
  //
  AliDebugClass(2,Form("SlatType %s Scale %e FlipX %d FlipY %d srcLine %d"
                       " destLine %d\n",slatType,scale,flipX,flipY,
                       srcLine,destLine));
  
  TString filename(AliMpFiles::SlatFilePath(kStationTrigger,slatType,
                                            planeType).Data());
  std::ifstream in(filename.Data());
  if (!in.good()) 
  {
    AliErrorClass(Form("Cannot read slat from %s",filename.Data()));
  }
  
  char line[80];
  
  while ( in.getline(line,80) )
  {
    TString sline(AliMpHelper::Normalize(line));

    if ( sline.Length() == 0 || sline[0] == '#' ) continue;
    
    Bool_t isKeywordThere = 
      sline.Contains(fgkKeywordPcb) || 
      sline.Contains(fgkKeywordLayer) ||
      sline.Contains(fgkKeywordScale) || 
      sline.Contains(fgkKeywordFlipX) || 
      sline.Contains(fgkKeywordFlipY);
    
    if ( !isKeywordThere ) 
    {
      AliErrorClass(Form("Got a line with no keyword : %s."
                         "That's not valid\n",line));
      continue; 
    }
    
    Double_t scale2;
    TString slatType2;
    
    Int_t isScaleLine = DecodeScaleLine(sline,scale2,slatType2);
    
    scale *= scale2;

    if ( isScaleLine < 0 )
    {
      AliFatalClass(Form("Syntax error near %s keyword\n",fgkKeywordScale.Data()));
    }
    else if ( isScaleLine > 0 && slatType2 != slatType )
    {
      ReadLines(slatType2.Data(),planeType,lines,scale,flipX,flipY,srcLine,destLine);
    }
    else    
    {
      Bool_t fx(kFALSE);
      Bool_t fy(kFALSE);
      Int_t isFlipLine = DecodeFlipLine(sline,slatType2,fx,fy);
      if ( isFlipLine )
      {
        if (fy)
        {
          srcLine = GetLine(slatType2);
          destLine = GetLine(slatType);
        }
        flipX |= fx;
        flipY |= fy;
        ReadLines(slatType2.Data(),planeType,lines,scale,flipX,flipY,srcLine,destLine);
      }
      else
      {
        lines.Add(new TObjString(sline.Data()));
      }
    }
  }
  
  in.close();
}
                                        
//_____________________________________________________________________________
void
AliMpTriggerReader::ReadLocalBoardMapping()
{
  // Reads the file that contains the mapping local board name <-> number

  TString filename(AliMpFiles::LocalTriggerBoardMapping());
  
  AliDebugClass(2,Form("Reading from %s\n",filename.Data()));

  fgLocalBoardMap.Delete();
  
  ifstream in(filename.Data());
  if (!in.good())
  {
    AliErrorClass(Form("Cannot read file %s\n",filename.Data()));    
  }
  else
  {
    char line[80];
    
    while ( in.getline(line,80) )
    {
      if ( line[0] == '#' ) continue;
      
      TString sline(line);
      if ( sline.Contains("Board") )
      {
        TObjArray* tokens = sline.Tokenize(' ');
        TString& number = ((TObjString*)(tokens->At(1)))->String();
        Int_t n = atoi(number.Data());
        if ( n == 0 ) continue;
        TString& name = ((TObjString*)(tokens->At(4)))->String();
        fgLocalBoardMap.Add(new TObjString(name), new TObjString(number));
        AliDebugClass(10,Form("Board %s has number %s\n",name.Data(),number.Data()));
        delete tokens;
      }
    }      
  }
  in.close();
}

//_____________________________________________________________________________
AliMpPCB*
AliMpTriggerReader::ReadPCB(const char* pcbType)
{ 
  //
  // Create a new AliMpPCB object, by reading it from file.
  //
  
  AliDebugClass(2,Form("pcbType=%s\n",pcbType));
  
  TString pcbName(pcbType);
  
  Ssiz_t pos = pcbName.First('x');

  Double_t scale = 1.0;
  
  if ( pos > 0 )
  {
    scale = TString(pcbName(pos+1,pcbName.Length()-pos-1)).Atof();
    pcbName = pcbName(0,pos);
  }
  
  std::ifstream in(AliMpFiles::SlatPCBFilePath(kStationTrigger,pcbName).Data());
  if (!in.good()) 
  {
    AliErrorClass(Form("Cannot open file for PCB %s",pcbName.Data()));
    return 0;
  }
 
  AliMpMotifReader reader(kStationTrigger,kNonBendingPlane); 
  // note that the nonbending
  // parameter is of no use for trigger, as far as reading motif is 
  // concerned, as all motifs are supposed to be in the same directory
  // (as they are shared by bending/non-bending planes).
     
  char line[80];
  
  const TString kSizeKeyword("SIZES");
  const TString kMotifKeyword("MOTIF");
  const TString kMotifSpecialKeyword("SPECIAL_MOTIF");
  
  AliMpPCB* pcb = 0;
  
  while ( in.getline(line,80) )
  {
    if ( line[0] == '#' ) continue;
    
    TString sline(line);
    
    if ( sline(0,kSizeKeyword.Length()) == kSizeKeyword )
    {
      std::istringstream sin(sline(kSizeKeyword.Length(),
                                   sline.Length()-kSizeKeyword.Length()-1).Data());
      float padSizeX = 0.0;
      float padSizeY = 0.0;
      float pcbSizeX = 0.0;
      float pcbSizeY = 0.0;
      sin >> padSizeX >> padSizeY >> pcbSizeX >> pcbSizeY;
      assert(pcb==0);
      pcb = new AliMpPCB(pcbType,padSizeX*scale,padSizeY*scale,
                         pcbSizeX*scale,pcbSizeY*scale);
    }
    
    if ( sline(0,kMotifSpecialKeyword.Length()) == kMotifSpecialKeyword )
    {
      std::istringstream sin(sline(kMotifSpecialKeyword.Length(),
                                   sline.Length()-kMotifSpecialKeyword.Length()).Data());
      TString sMotifSpecial;
      TString sMotifType;
      sin >> sMotifSpecial >> sMotifType;
      
      AliMpMotifType* motifType = reader.BuildMotifType(sMotifType);
      AliMpMotifSpecial* specialMotif = 
        reader.BuildMotifSpecial(sMotifSpecial,motifType,scale);
      
      assert(pcb==0);      
      pcb = new AliMpPCB(pcbType,specialMotif);
    }
    
    if ( sline(0,kMotifKeyword.Length()) == kMotifKeyword )
    {
      std::istringstream sin(sline(kMotifKeyword.Length(),
                                   sline.Length()-kMotifKeyword.Length()).Data());
      TString sMotifType;
      int ix;
      int iy;
      sin >> sMotifType >> ix >> iy;
      
      AliMpMotifType* motifType = reader.BuildMotifType(sMotifType.Data());
      
      assert(pcb!=0);
      pcb->Add(motifType,ix,iy);
    }
  }
  
  in.close();
  
  fgPCBMap.Add(new TObjString(pcbType),pcb);
  return pcb;
}

//_____________________________________________________________________________
AliMpTrigger*
AliMpTriggerReader::ReadSlat(const char* slatType, AliMpPlaneType planeType)
{
  //
  // Create a new AliMpTrigger object, by reading it from file.
  //

  Double_t scale = 1.0;
  Bool_t flipX = kFALSE;
  Bool_t flipY = kFALSE;
  TList lines;
  Int_t srcLine(-1);
  Int_t destLine(-1);
  
  // Read the file and its include (if any) and store the result
  // in a TObjArray of TObjStrings.
  ReadLines(slatType,planeType,lines,scale,flipX,flipY,srcLine,destLine);

  // Here some more sanity checks could be done.
  // For the moment we only insure that the first line contains 
  // a layer keyword.
  TString& firstLine = ((TObjString*)lines.First())->String();
  if ( !IsLayerLine(firstLine) ) 
  {
    std::ostringstream s;
    s << fgkKeywordLayer;
    lines.AddFirst(new TObjString(s.str().c_str()));
  }
  
  AliDebugClass(2,Form("Scale=%g\n",scale));
  
  FlipLines(lines,flipX,flipY,srcLine,destLine);
  
  // Now splits the lines in packets corresponding to different layers 
  // (if any), and create sub-slats.
  TObjArray layers;
  Int_t ilayer(-1);
  TIter it(&lines);
  TObjString* osline;
  
  while ( ( osline = (TObjString*)it.Next() ) )
  {
    TString& s = osline->String();
    if ( IsLayerLine(s) )
    {
      layers.Add(new TList);
      ++ilayer;
    }
    else
    {
      ((TList*)layers.At(ilayer))->Add(new TObjString(s));
    }
  }

  AliDebugClass(2,Form("nlayers=%d\n",layers.GetEntriesFast()));

  AliMpTrigger* triggerSlat = new AliMpTrigger(slatType, planeType);
    
  for ( Int_t ilayer = 0; ilayer < layers.GetEntriesFast(); ++ilayer )
  {
    TList& lines = *((TList*)layers.At(ilayer));
    std::ostringstream slatName;
    slatName << slatType << "-LAYER" << ilayer;
    AliMpSlat* slat = BuildSlat(slatName.str().c_str(),planeType,lines,scale);
    if ( slat )
    {
      triggerSlat->AdoptLayer(slat);
    }
    else
    {
      AliErrorClass(Form("Could not read %s\n",slatName.str().c_str()));
      delete triggerSlat;
      return 0;
    }
  }
  
  layers.SetOwner(kTRUE);
  layers.Delete();
  
  return triggerSlat;
}

//_____________________________________________________________________________
void
AliMpTriggerReader::Reset()
{
  // Resets the PCB internal map
  fgPCBMap.Delete();
}
