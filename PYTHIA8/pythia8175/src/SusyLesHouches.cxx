// SusyLesHouches.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// Main authors of this file: N. Desai, P. Skands
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#include "SusyLesHouches.h"

// GZIP support.
#ifdef GZIPSUPPORT

// For GCC versions >= 4.6, can switch off shadow warnings.
#if (defined GZIPSUPPORT && ((__GNUC__ * 100) + __GNUC_MINOR__) >= 406)
#pragma GCC diagnostic ignored "-Wshadow"
#endif

// Boost includes.
#include "boost/iostreams/filtering_stream.hpp"
#include "boost/iostreams/filter/gzip.hpp"

// Switch shadow warnings back on.
#if (defined GZIPSUPPORT && ((__GNUC__ * 100) + __GNUC_MINOR__) >= 406)
#pragma GCC diagnostic warning "-Wshadow"
#endif

#endif // GZIPSUPPORT

namespace Pythia8 {

//==========================================================================

// The SusyLesHouches class.

//==========================================================================

// Main routine to read in SLHA and LHEF+SLHA files

int SusyLesHouches::readFile(string slhaFileIn, int verboseIn, 
  bool useDecayIn) {

  // Copy inputs to local
  slhaFile = slhaFileIn;
  verbose  = verboseIn;
  useDecay = useDecayIn;

  // Check that input file is OK.
  int iFailFile=0;
  const char* cstring = slhaFile.c_str();

// Construct istream without gzip support.
#ifndef GZIPSUPPORT
  ifstream file(cstring);

// Construct istream with gzip support.
#else
  boost::iostreams::filtering_istream file;
  ifstream fileBase(cstring);

  // Pass along the 'good()' flag, so code elsewhere works unmodified.
  if (!fileBase.good()) file.setstate(ios_base::badbit);

  // Check filename ending to decide which filters to apply.
  else {
    const char *last = strrchr(cstring, '.');
    if (last && strncmp(last, ".gz", 3) == 0)
      file.push(boost::iostreams::gzip_decompressor());
    file.push(fileBase);
  }
#endif

  // Exit if input file not found. Else print file name.
  if (!file.good()) {
    message(2,"readFile",slhaFile+" not found",0);
    return -1;
    slhaRead=false;
  }  
  if (verbose >= 3) {
    message(0,"readFile","parsing "+slhaFile,0);
    filePrinted = true;
  }

  // Array of particles read in.
  vector<int> idRead;

  //Initial values for read-in variables.  
  slhaRead=true;
  lhefRead=false;
  lhefSlha=false;
  bool foundSlhaTag = false;
  bool xmlComment   = false;
  bool decayPrinted = false;
  string line="";
  string blockIn="";
  string decay="";
  string comment="";
  string blockName="";
  string nameNow="";
  int idNow=0;
  double width=0.0;

  //Initialize line counter
  int iLine=0;

  // Read in one line at a time.
  while ( getline(file, line) ) {
    iLine++;

    //Rewrite string in lowercase
    for (unsigned int i=0;i<line.length();i++) line[i]=tolower(line[i]);

    // Remove extra blanks 
    while (line.find("  ") != string::npos) line.erase( line.find("  "), 1);

    //Detect whether read-in is from a Les Houches Event File (LHEF).
    if (line.find("<leshouches") != string::npos 
        || line.find("<slha") != string::npos) {
      lhefRead=true;
    }

    // If LHEF
    if (lhefRead) {      
      //Ignore XML comments (only works for whole lines so far)
      if (line.find("-->") != string::npos) {
        xmlComment = false;
      }
      else if (xmlComment) continue;
      else if (line.find("<!--") != string::npos) {
        xmlComment = true;
      } 
      //Detect when <slha> tag reached.
      if (line.find("<slha") != string::npos) {
        lhefSlha     = true;
        foundSlhaTag = true;
        //Print header if not already done
        if (! headerPrinted) printHeader();
      }
      //Stop looking when </header> or <init> tag reached
      if (line.find("</header>") != string::npos ||
          line.find("<init") != string::npos) {      
        if (!foundSlhaTag) return 101;
        break;
      }
      //If <slha> tag not yet reached, skip
      if (!lhefSlha) continue;
    }

    //Ignore comment lines with # as first character
    if (line.find("#") == 0) continue;

    //Ignore empty lines
    if (line.size() == 0) continue;
    if (line.size() == 1 && line.substr(0,1) == " ") continue;

    //Move comment to separate string
    if (line.find("#") != string::npos) {
      if (line.find("#") + 1 < line.length() )
	comment = line.substr(line.find("#")+1,line.length()-line.find("#")-2);
      else 
	comment = "";
      line.erase(line.find("#"),line.length()-line.find("#")-1);
    }

    // Remove blanks before and after an = sign.
    while (line.find(" =") != string::npos) line.erase( line.find(" ="), 1);
    while (line.find("= ") != string::npos) line.erase( line.find("= ")+1, 1);

    //New block. 
    if (line.find("block") <= 1) { 

      //Print header if not already done
      if (! headerPrinted) printHeader();

      blockIn=line ;       
      decay="";
      int nameBegin=6 ;
      int nameEnd=blockIn.find(" ",7);
      blockName=blockIn.substr(nameBegin,nameEnd-nameBegin);
      
      // Copy input file as generic blocks (containing strings)
      // (more will be done with SLHA1 & 2 specific blocks below, this is 
      //  just to make sure we have a complete copy of the input file, 
      //  including also any unknown/user/generic blocks)
      LHgenericBlock gBlock;
      genericBlocks[blockName]=gBlock;

      // QNUMBERS blocks (cf. arXiv:0712.3311 [hep-ph])
      if (blockIn.find("qnumbers") != string::npos) {
	// Extract ID code for new particle
	int pdgBegin=blockIn.find(" ",7)+1;
	int pdgEnd=blockIn.find(" ",pdgBegin);
	string pdgString = blockIn.substr(pdgBegin,pdgEnd-pdgBegin);
	istringstream linestream(pdgString);
	// Create and add new block with this code as zero'th entry
	LHblock<int> newQnumbers;
	newQnumbers.set(0,linestream);
	qnumbers.push_back(newQnumbers);	
	// Default name: PDG code
	string defName, defAntiName, newName, newAntiName;	
	ostringstream idStream;
	idStream<<newQnumbers(0);
	defName     = idStream.str();
	defAntiName = "-"+defName;
	newName     = defName;
	newAntiName = defAntiName;
	// Attempt to extract names from comment string
	if (comment.length() >= 1) {
	  int firstCommentBeg(0), firstCommentEnd(0);
	  if ( comment.find(" ") == 0) firstCommentBeg = 1;	  
	  if ( comment.find(" ",firstCommentBeg+1) == string::npos)
	    firstCommentEnd = comment.length();
	  else 
	    firstCommentEnd = comment.find(" ",firstCommentBeg+1);
	  if (firstCommentEnd > firstCommentBeg) 
	    newName = comment.substr(firstCommentBeg,
				     firstCommentEnd-firstCommentBeg);
	  // Now see if there is a separate name for antiparticle
	  int secondCommentBeg(firstCommentEnd+1), secondCommentEnd(0);
	  if (secondCommentBeg < int(comment.length())) { 
	    if ( comment.find(" ",secondCommentBeg+1) == string::npos)
	      secondCommentEnd = comment.length();
	    else 
	      secondCommentEnd = comment.find(" ",secondCommentBeg+1);	    
	    if (secondCommentEnd > secondCommentBeg) 
	      newAntiName = comment.substr(secondCommentBeg,
					   secondCommentEnd-secondCommentBeg);
	  }	  
	} 
	// If name given without specific antiname, set antiname to ""
	if (newName != defName && newAntiName == defAntiName) newAntiName = "";
	qnumbersName.push_back(newName);
	qnumbersAntiName.push_back(newAntiName);
	if (pdgString != newName) {
	  message(0,"readFile","storing QNUMBERS for id = "+pdgString+" "
		  +newName+" "+newAntiName,iLine);
	} else {
	  message(0,"readFile","storing QNUMBERS for id = "+pdgString,iLine);
	}
      }

      //Find Q=... for DRbar running blocks
      if (blockIn.find("q=") != string::npos) {
        int qbegin=blockIn.find("q=")+2;
        istringstream qstream(blockIn.substr(qbegin,blockIn.length()));
        double q=0.0;
        qstream >> q;
        if (qstream) {
          // SLHA1 running blocks
          if (blockName=="hmix") hmix.setq(q);
          if (blockName=="yu") yu.setq(q);
          if (blockName=="yd") yd.setq(q);
          if (blockName=="ye") ye.setq(q);
          if (blockName=="au") au.setq(q);
          if (blockName=="ad") ad.setq(q);
          if (blockName=="ae") ae.setq(q);
          if (blockName=="msoft") msoft.setq(q);
          if (blockName=="gauge") gauge.setq(q);
          // SLHA2 running blocks
          if (blockName=="vckm") vckm.setq(q);
          if (blockName=="upmns") upmns.setq(q);
          if (blockName=="msq2") msq2.setq(q);
          if (blockName=="msu2") msu2.setq(q);
          if (blockName=="msd2") msd2.setq(q);
          if (blockName=="msl2") msl2.setq(q);
          if (blockName=="mse2") mse2.setq(q);
          if (blockName=="tu") tu.setq(q);
          if (blockName=="td") td.setq(q);
          if (blockName=="te") te.setq(q);
          if (blockName=="rvlamlle") rvlamlle.setq(q);
          if (blockName=="rvlamlqd") rvlamlqd.setq(q);
          if (blockName=="rvlamudd") rvlamudd.setq(q);
          if (blockName=="rvtlle") rvtlle.setq(q);
          if (blockName=="rvtlqd") rvtlqd.setq(q);
          if (blockName=="rvtudd") rvtudd.setq(q);
          if (blockName=="rvkappa") rvkappa.setq(q);
          if (blockName=="rvd") rvd.setq(q);
          if (blockName=="rvm2lh1") rvm2lh1.setq(q);
          if (blockName=="rvsnvev") rvsnvev.setq(q);          
          if (blockName=="imau") imau.setq(q);
          if (blockName=="imad") imad.setq(q);
          if (blockName=="imae") imae.setq(q);
          if (blockName=="imhmix") imhmix.setq(q);
          if (blockName=="immsoft") immsoft.setq(q);
          if (blockName=="imtu") imtu.setq(q);
          if (blockName=="imtd") imtd.setq(q);
          if (blockName=="imte") imte.setq(q);
          if (blockName=="imvckm") imvckm.setq(q);
          if (blockName=="imupmns") imupmns.setq(q);
          if (blockName=="immsq2") immsq2.setq(q);
          if (blockName=="immsu2") immsu2.setq(q);
          if (blockName=="immsd2") immsd2.setq(q);
          if (blockName=="immsl2") immsl2.setq(q);
          if (blockName=="immse2") immse2.setq(q);
        };
      };
      
      //Skip to next line.
      continue ; 
      
    } 

    //New decay table
    else if (line.find("decay") <= 1) {

      // Print header if not already done
      if (! headerPrinted) printHeader();

      // If previous had zero length, print now
      if (decay != "" && ! decayPrinted) {
        if (verbose >= 2) message(0,"readFile","reading  WIDTH for "+nameNow
                +" (but no decay channels found)",0);
      }

      //Set decay block name
      decay=line;
      blockIn="";
      int nameBegin=6 ;
      int nameEnd=decay.find(" ",7);
      nameNow=decay.substr(nameBegin,nameEnd-nameBegin);
      
      //Extract PDG code and width
      istringstream dstream(nameNow);
      dstream >> idNow;

      //Ignore decay if decay table read-in switched off
      if( !useDecay ) {
	decay = "";
	message(0,"readFile","ignoring DECAY table for "+nameNow
		+" (DECAY read-in switched off)",iLine);
	continue;
      }

      if (dstream) {
        string widthName=decay.substr(nameEnd+1,decay.length());
        istringstream wstream(widthName);
        wstream >> width;
        if (wstream) {
          // Set 
          decays.push_back(LHdecayTable(idNow,width));          
          decayIndices[idNow]=decays.size()-1;
          //Set PDG code and width
          if (width <= 0.0) {
            string endComment="";
            if (width < -1e-6) {
              endComment="(forced width < 0 to zero)";
            }
            if (verbose >= 2)
              message(0,"readFile","reading  stable particle "+nameNow
                      +" "+endComment,0);
            width=0.0;
            decayPrinted = true;
            decays[decayIndices[idNow]].setWidth(width);
          } else {
            decayPrinted = false;
          }
        } else {
          if (verbose >= 2) 
            message(0,"readFile","ignoring DECAY table for "+nameNow
                    +" (read failed)",iLine);
          decayPrinted = true;
          width=0.0;
          decay="";
          continue;
        }
      }
      else {
        message(0,"readFile",
                    "PDG Code unreadable. Ignoring this DECAY block",iLine);
        decayPrinted = true;
        decay="";
        continue;
      }

      //Skip to next line
      continue ;
    }

    //Switch off SLHA read-in via LHEF if outside <slha> tag.
    else if (line.find("</slha>") != string::npos) {
      lhefSlha=false;
      blockIn="";
      decay="";
      continue;
    }

    //Skip not currently reading block data lines.
    if (blockIn != "") {

      // Replace an equal sign by a blank to make parsing simpler.
      while (line.find("=") != string::npos) {
        int firstEqual = line.find_first_of("=");
        line.replace(firstEqual, 1, " ");   
      };
    
      //Parse data lines within given block
      //Constructed explicitly so that each block can have its own types and
      //own rules defined. For extra user blocks, just add more recognized 
      //blockNames at the end and implement user defined rules accordingly.
      //string comment = line.substr(line.find("#"),line.length());    
      int ifail=-2;
      istringstream linestream(line);

      // Read line in QNUMBERS block, add to end of qnumbers vector
      if (blockName == "qnumbers") {
	int iEnd = qnumbers.size()-1;
	if (iEnd >= 0) ifail = qnumbers[iEnd].set(linestream);
	else ifail = -1;
      }

      // MODEL
      else if (blockName == "modsel") {
        int i;
        linestream >> i; 
        if (linestream) {
          if (i == 12) {ifail=modsel12.set(0,linestream);} 
          else if (i == 21) {ifail=modsel21.set(0,linestream);}
          else {ifail=modsel.set(i,linestream);};}
        else {
          ifail = -1;}
      };
      if (blockName == "minpar") ifail=minpar.set(linestream); 
      if (blockName == "sminputs") ifail=sminputs.set(linestream);
      if (blockName == "extpar") ifail=extpar.set(linestream);
      if (blockName == "qextpar") ifail=qextpar.set(linestream);
      //FLV
      if (blockName == "vckmin") ifail=vckmin.set(linestream);
      if (blockName == "upmnsin") ifail=upmnsin.set(linestream);
      if (blockName == "msq2in") ifail=msq2in.set(linestream);
      if (blockName == "msu2in") ifail=msu2in.set(linestream);
      if (blockName == "msd2in") ifail=msd2in.set(linestream);
      if (blockName == "msl2in") ifail=msl2in.set(linestream);
      if (blockName == "mse2in") ifail=mse2in.set(linestream);
      if (blockName == "tuin") ifail=tuin.set(linestream);
      if (blockName == "tdin") ifail=tdin.set(linestream);
      if (blockName == "tein") ifail=tein.set(linestream);
      //RPV
      if (blockName == "rvlamllein") ifail=rvlamllein.set(linestream);
      if (blockName == "rvlamlqdin") ifail=rvlamlqdin.set(linestream);
      if (blockName == "rvlamuddin") ifail=rvlamuddin.set(linestream);
      if (blockName == "rvtllein") ifail=rvtllein.set(linestream);
      if (blockName == "rvtlqdin") ifail=rvtlqdin.set(linestream);
      if (blockName == "rvtuddin") ifail=rvtuddin.set(linestream);
      if (blockName == "rvkappain") ifail=rvkappain.set(linestream);
      if (blockName == "rvdin") ifail=rvdin.set(linestream);
      if (blockName == "rvm2lh1in") ifail=rvm2lh1in.set(linestream);
      if (blockName == "rvsnvevin") ifail=rvsnvevin.set(linestream);
      //CPV 
      if (blockName == "imminpar") ifail=imminpar.set(linestream);
      if (blockName == "imextpar") ifail=imextpar.set(linestream);
      //CPV +FLV
      if (blockName == "immsq2in") ifail=immsq2in.set(linestream);
      if (blockName == "immsu2in") ifail=immsu2in.set(linestream);
      if (blockName == "immsd2in") ifail=immsd2in.set(linestream);
      if (blockName == "immsl2in") ifail=immsl2in.set(linestream);
      if (blockName == "immse2in") ifail=immse2in.set(linestream);
      if (blockName == "imtuin") ifail=imtuin.set(linestream);
      if (blockName == "imtdin") ifail=imtdin.set(linestream);
      if (blockName == "imtein") ifail=imtein.set(linestream);
      //Info:
      if (blockName == "spinfo" || blockName=="dcinfo") {
        int i;
        string entry;
        linestream >> i >> entry;
        string blockStr="RGE";
        if (blockName=="dcinfo") blockStr="DCY";

        if (linestream) {
          if ( i == 3 ) {
            string warning=line.substr(line.find("3")+1,line.length());
            message(1,"readFile","(from "+blockStr+" program): "+warning,0);
            if (blockName == "spinfo") spinfo3.set(warning);
            else dcinfo3.set(warning);
          } else if ( i == 4 ) {
            string error=line.substr(line.find("4")+1,line.length());
            message(2,"readFile","(from "+blockStr+" program): "+error,0);
            if (blockName == "spinfo") spinfo4.set(error);
            else dcinfo4.set(error);
          } else {
            //Rewrite string in uppercase
            for (unsigned int j=0;j<entry.length();j++) 
              entry[j]=toupper(entry[j]);
            ifail=(blockName=="spinfo") ? spinfo.set(i,entry)
              : dcinfo.set(i,entry);
          };
        } else {
          ifail=-1;
        };
      };
      //SPECTRUM
      //Pole masses
      if (blockName == "mass") ifail=mass.set(linestream);

      //Mixing
      if (blockName == "alpha") ifail=alpha.set(linestream,false);
      if (blockName == "stopmix") ifail=stopmix.set(linestream);
      if (blockName == "sbotmix") ifail=sbotmix.set(linestream);
      if (blockName == "staumix") ifail=staumix.set(linestream);
      if (blockName == "nmix") ifail=nmix.set(linestream);
      if (blockName == "umix") ifail=umix.set(linestream);
      if (blockName == "vmix") ifail=vmix.set(linestream);
      //FLV
      if (blockName == "usqmix") ifail=usqmix.set(linestream);
      if (blockName == "dsqmix") ifail=dsqmix.set(linestream);
      if (blockName == "selmix") ifail=selmix.set(linestream);
      if (blockName == "snumix") ifail=snumix.set(linestream);
      if (blockName == "snsmix") ifail=snsmix.set(linestream);
      if (blockName == "snamix") ifail=snamix.set(linestream);
      //RPV
      if (blockName == "rvnmix") ifail=rvnmix.set(linestream);
      if (blockName == "rvumix") ifail=rvumix.set(linestream);
      if (blockName == "rvvmix") ifail=rvvmix.set(linestream);
      if (blockName == "rvhmix") ifail=rvhmix.set(linestream);
      if (blockName == "rvamix") ifail=rvamix.set(linestream);
      if (blockName == "rvlmix") ifail=rvlmix.set(linestream);
      //CPV
      if (blockName == "cvhmix") ifail=cvhmix.set(linestream);
      if (blockName == "imcvhmix") ifail=imcvhmix.set(linestream);
      //CPV + FLV
      if (blockName == "imusqmix") ifail=imusqmix.set(linestream);
      if (blockName == "imdsqmix") ifail=imdsqmix.set(linestream);
      if (blockName == "imselmix") ifail=imselmix.set(linestream);
      if (blockName == "imsnumix") ifail=imsnumix.set(linestream);
      if (blockName == "imnmix") ifail=imnmix.set(linestream);
      if (blockName == "imumix") ifail=imumix.set(linestream);
      if (blockName == "imvmix") ifail=imvmix.set(linestream);
      //NMSSM
      if (blockName == "nmhmix") ifail=nmhmix.set(linestream);
      if (blockName == "nmamix") ifail=nmamix.set(linestream);
      if (blockName == "nmnmix") ifail=nmnmix.set(linestream);
      
      //DRbar Lagrangian parameters
      if (blockName == "gauge") ifail=gauge.set(linestream);      
      if (blockName == "yu") ifail=yu.set(linestream);
      if (blockName == "yd") ifail=yd.set(linestream);
      if (blockName == "ye") ifail=ye.set(linestream);
      if (blockName == "au") ifail=au.set(linestream);
      if (blockName == "ad") ifail=ad.set(linestream);
      if (blockName == "ae") ifail=ae.set(linestream);
      if (blockName == "hmix") ifail=hmix.set(linestream);
      if (blockName == "msoft") ifail=msoft.set(linestream);
      //FLV
      if (blockName == "vckm") ifail=vckm.set(linestream);
      if (blockName == "upmns") ifail=upmns.set(linestream);
      if (blockName == "msq2") ifail=msq2.set(linestream);
      if (blockName == "msu2") ifail=msu2.set(linestream);
      if (blockName == "msd2") ifail=msd2.set(linestream);
      if (blockName == "msl2") ifail=msl2.set(linestream);
      if (blockName == "mse2") ifail=mse2.set(linestream);
      if (blockName == "tu") ifail=tu.set(linestream);
      if (blockName == "td") ifail=td.set(linestream);
      if (blockName == "te") ifail=te.set(linestream);
      //RPV
      if (blockName == "rvlamlle") ifail=rvlamlle.set(linestream);
      if (blockName == "rvlamlqd") ifail=rvlamlqd.set(linestream);
      if (blockName == "rvlamudd") ifail=rvlamudd.set(linestream);
      if (blockName == "rvtlle") ifail=rvtlle.set(linestream);
      if (blockName == "rvtlqd") ifail=rvtlqd.set(linestream);
      if (blockName == "rvtudd") ifail=rvtudd.set(linestream);
      if (blockName == "rvkappa") ifail=rvkappa.set(linestream);
      if (blockName == "rvd") ifail=rvd.set(linestream);
      if (blockName == "rvm2lh1") ifail=rvm2lh1.set(linestream);
      if (blockName == "rvsnvev") ifail=rvsnvev.set(linestream);
      //CPV
      if (blockName == "imau") ifail=imau.set(linestream);
      if (blockName == "imad") ifail=imad.set(linestream);
      if (blockName == "imae") ifail=imae.set(linestream);
      if (blockName == "imhmix") ifail=imhmix.set(linestream);
      if (blockName == "immsoft") ifail=immsoft.set(linestream);
      //CPV+FLV
      if (blockName == "imvckm") ifail=imvckm.set(linestream);
      if (blockName == "imupmns") ifail=imupmns.set(linestream);
      if (blockName == "immsq2") ifail=immsq2.set(linestream);
      if (blockName == "immsu2") ifail=immsu2.set(linestream);
      if (blockName == "immsd2") ifail=immsd2.set(linestream);
      if (blockName == "immsl2") ifail=immsl2.set(linestream);
      if (blockName == "immse2") ifail=immse2.set(linestream);
      if (blockName == "imtu") ifail=imtu.set(linestream);
      if (blockName == "imtd") ifail=imtd.set(linestream);
      if (blockName == "imte") ifail=imte.set(linestream);
      //NMSSM
      if (blockName == "nmssmrun") ifail=nmssmrun.set(linestream);      

      //Diagnostics
      if (ifail != 0) { 
        if (ifail == -2 && !genericBlocks[blockName].exists() ) {
          message(0,"readFile","storing non-SLHA(2) block: "+blockName,iLine);
        };
        if (ifail == -1) {
          message(1,"readFile","read error or empty line",iLine);        
        };
        if (ifail == 1) {
          message(0,"readFile",blockName+" existing entry overwritten",iLine);
        };
      }

      // Add line to generic block (carbon copy of input structure)
      // NB: do not save empty lines, defined as having length <= 1
      if (line.size() >= 2) {
	genericBlocks[blockName].set(line);
      }
	
    } 

    // Decay table read-in
    else if (decay != "") {
      if (! decayPrinted) {
        if (verbose >= 2) 
          message(0,"readFile","reading  DECAY table for "+nameNow,0);
        decayPrinted = true;
      }
      double brat;
      bool ok=true;
      int nDa = 0;
      vector<int> idDa;
      istringstream linestream(line);
      linestream >> brat;
      if (! linestream) ok = false;
      if (ok) linestream >> nDa;
      if (! linestream) ok = false;
      else {
        for (int i=0; i<nDa; i++) {
          int idThis;
          linestream >> idThis;
          if (! linestream) {
            ok = false;
            break;
          }
          idDa.push_back(idThis);
        }
      }

      // Stop reading decay channels if not consistent.
      if (!ok || nDa < 2) {
        message(1,"readFile","read error or empty line",iLine);
         
      // Append decay channel.
      } else {
        decays[decayIndices[idNow]].addChannel(brat,nDa,idDa);
      }
    }
  };

  //Print footer 
  printFooter();

  //Return 0 if read-in successful 
  if ( lhefRead && !foundSlhaTag) { 
    return 102; 
  }
  else return iFailFile;
    
}

//--------------------------------------------------------------------------

// Print a header with information on version, last date of change, etc.

void SusyLesHouches::printHeader() {
  if (verbose == 0) return;
  setprecision(3);
  if (! headerPrinted) {
    cout << " *-----------------------  SusyLesHouches SUSY/BSM"
         << " Interface  ------------------------*\n";
    message(0,"","Last Change 01 Aug 2012 - P. Skands",0);
    if (!filePrinted) {
      message(0,"","Parsing: "+slhaFile,0);
      filePrinted=true;
    }
    headerPrinted=true;
  }
}

//--------------------------------------------------------------------------

// Print a footer

void SusyLesHouches::printFooter() {
  if (verbose == 0) return;
  if (! footerPrinted) {
    //    cout << " *"<<endl;
    cout << " *-----------------------------------------------------"
         << "-------------------------------*\n";
    footerPrinted=true;
    //    headerPrinted=false; 
  }
}

//--------------------------------------------------------------------------

// Print the current spectrum on stdout.
// Not yet fully implemented.

void SusyLesHouches::printSpectrum(int ifail) {

  // Exit if output switched off
  if (verbose <= 0) return;

  // Print header if not already done
  if (! headerPrinted) printHeader();
  message(0,"","");

  // Print Calculator and File name
  if (slhaRead) {
    message(0,"","  Spectrum Calculator was:   "+spinfo(1)+"   version: "
      +spinfo(2));
    if (lhefRead) message(0,"","  Read <slha> spectrum from: "+slhaFile);
    else message(0,"","  Read SLHA spectrum from: "+slhaFile);
  }

  // Failed?
  if (ifail < 0) {
    message(0,"","  Check revealed problems. Only using masses.");
  }

  // gluino
  message(0,"","");
  cout<<" |  ~g                  m"<<endl; 
  cout<<setprecision(3)<<" |     1000021 "<<setw(10)<<
      ( (mass(2000003) > 1e7) ? scientific : fixed)<<mass(1000021)<<endl;

  // d squarks
  message(0,"","");
  cout << " |  ~d                  m     ~dL     ~sL     ~bL"
       << "     ~dR     ~sR     ~bR" << endl;

  cout<<setprecision(3) <<" |     1000001 "<<setw(10)<<
    ( (mass(1000001) > 1e7) ? scientific : fixed)<<mass(1000001)<<fixed<<"  ";
  for (int icur=1;icur<=6;icur++) cout<<setw(6)<<dsqmix(1,icur)<<"  ";

  cout<<endl<<" |     1000003 "<<setw(10)<<
    ( (mass(1000003) > 1e7) ? scientific : fixed)<<mass(1000003)<<fixed<<"  ";
  for (int icur=1;icur<=6;icur++) cout<<setw(6)<<dsqmix(2,icur)<<"  ";

  cout<<endl<<" |     1000005 "<<setw(10)<<
    ( (mass(1000005) > 1e7) ? scientific : fixed)<<mass(1000005)<<fixed<<"  ";
  for (int icur=1;icur<=6;icur++) cout<<setw(6)<<dsqmix(3,icur)<<"  ";

  cout<<endl<<" |     2000001 "<<setw(10)<<
    ( (mass(2000001) > 1e7) ? scientific : fixed)<<mass(2000001)<<fixed<<"  ";
  for (int icur=1;icur<=6;icur++) cout<<setw(6)<<dsqmix(4,icur)<<"  ";

  cout<<endl<<" |     2000003 "<<setw(10)<<
    ( (mass(2000003) > 1e7) ? scientific : fixed)<<mass(2000003)<<fixed<<"  ";
  for (int icur=1;icur<=6;icur++) cout<<setw(6)<<dsqmix(5,icur)<<"  ";

  cout<<endl<<" |     2000005 "<<setw(10)<<
    ( (mass(2000005) > 1e7) ? scientific : fixed)<<mass(2000005)<<fixed<<"  ";
  for (int icur=1;icur<=6;icur++) cout<<setw(6)<<dsqmix(6,icur)<<"  ";

  cout<<endl;
  
  // u squarks
  message(0,"","");
  cout << " |  ~u                  m     ~uL     ~cL     ~tL"
       << "     ~uR     ~cR     ~tR" << endl; 

  cout<<setprecision(3)<<" |     1000002 "<<setw(10)<<
    ( (mass(1000002) > 1e7) ? scientific : fixed)<<mass(1000002)<<fixed<<"  ";
  for (int icur=1;icur<=6;icur++) cout<<setw(6)<<usqmix(1,icur)<<"  ";

  cout<<endl<<" |     1000004 "<<setw(10)<<
    ( (mass(1000004) > 1e7) ? scientific : fixed)<<mass(1000004)<<fixed<<"  ";
  for (int icur=1;icur<=6;icur++) cout<<setw(6)<<usqmix(2,icur)<<"  ";

  cout<<endl<<" |     1000006 "<<setw(10)<<
    ( (mass(1000006) > 1e7) ? scientific : fixed)<<mass(1000006)<<fixed<<"  "; 
  for (int icur=1;icur<=6;icur++) cout<<setw(6)<<usqmix(3,icur)<<"  ";

  cout<<endl<<" |     2000002 "<<setw(10)<<
    ( (mass(2000002) > 1e7) ? scientific : fixed)<<mass(2000002)<<fixed<<"  "; 
  for (int icur=1;icur<=6;icur++) cout<<setw(6)<<usqmix(4,icur)<<"  ";

  cout<<endl<<" |     2000004 "<<setw(10)<<
    ( (mass(2000004) > 1e7) ? scientific : fixed)<<mass(2000004)<<fixed<<"  " ;
  for (int icur=1;icur<=6;icur++) cout<<setw(6)<<usqmix(5,icur)<<"  ";

  cout<<endl<<" |     2000006 "<<setw(10)<<
    ( (mass(2000006) > 1e7) ? scientific : fixed)<<mass(2000006)<<fixed<<"  ";
  for (int icur=1;icur<=6;icur++) cout<<setw(6)<<usqmix(6,icur)<<"  ";

  cout<<endl;

  // Charged scalars (sleptons)
  message(0,"","");

  // R-conserving:
  if (modsel(4) < 1) {
    cout << " |  ~e                  m     ~eL    ~muL   ~tauL"
	 << "     ~eR    ~muR   ~tauR" << endl; 

    cout<<setprecision(3)<<" |     1000011 "<<setw(10)<<
      ( (mass(1000011) > 1e7) ? scientific : fixed)<<mass(1000011)<<fixed<<"  ";
    for (int icur=1;icur<=6;icur++) cout<<setw(6)<<selmix(1,icur)<<"  ";
    
    cout<<endl<<" |     1000013 "<<setw(10)<<
      ( (mass(1000013) > 1e7) ? scientific : fixed)<<mass(1000013)<<fixed<<"  ";
    for (int icur=1;icur<=6;icur++) cout<<setw(6)<<selmix(2,icur)<<"  ";
    
    cout<<endl<<" |     1000015 "<<setw(10)<<
      ( (mass(1000015) > 1e7) ? scientific : fixed)<<mass(1000015)<<fixed<<"  "; 
    for (int icur=1;icur<=6;icur++) cout<<setw(6)<<selmix(3,icur)<<"  ";
    
    cout<<endl<<" |     2000011 "<<setw(10)<<
      ( (mass(2000011) > 1e7) ? scientific : fixed)<<mass(2000011)<<fixed<<"  "; 
    for (int icur=1;icur<=6;icur++) cout<<setw(6)<<selmix(4,icur)<<"  ";
    
    cout<<endl<<" |     2000013 "<<setw(10)<<
      ( (mass(2000013) > 1e7) ? scientific : fixed)<<mass(2000013)<<fixed<<"  " ;
    for (int icur=1;icur<=6;icur++) cout<<setw(6)<<selmix(5,icur)<<"  ";
    
    cout<<endl<<" |     2000015 "<<setw(10)<<
      ( (mass(2000015) > 1e7) ? scientific : fixed)<<mass(2000015)<<fixed<<"  ";
    for (int icur=1;icur<=6;icur++) cout<<setw(6)<<selmix(6,icur)<<"  ";
  }

  // R-violating
  else {
    cout << " |  H-/~e               m     H1-     H2-     ~eL    ~muL   ~tauL"
	 << "     ~eR    ~muR   ~tauR" << endl; 

    cout<<setprecision(3)<<" |         -37 "<<setw(10)<<
      ( (mass(37) > 1e7) ? scientific : fixed)<<mass(37)<<fixed<<"  ";
    for (int icur=1;icur<=8;icur++) cout<<setw(6)<<rvlmix(1,icur)<<"  ";
    
    cout<<endl<<" |     1000011 "<<setw(10)<<
      ( (mass(1000011) > 1e7) ? scientific : fixed)<<mass(1000011)<<fixed<<"  ";
    for (int icur=1;icur<=8;icur++) cout<<setw(6)<<rvlmix(2,icur)<<"  ";
    
    cout<<endl<<" |     1000013 "<<setw(10)<<
      ( (mass(1000013) > 1e7) ? scientific : fixed)<<mass(1000013)<<fixed<<"  ";
    for (int icur=1;icur<=8;icur++) cout<<setw(6)<<rvlmix(3,icur)<<"  ";
    
    cout<<endl<<" |     1000015 "<<setw(10)<<
      ( (mass(1000015) > 1e7) ? scientific : fixed)<<mass(1000015)<<fixed<<"  "; 
    for (int icur=1;icur<=8;icur++) cout<<setw(6)<<rvlmix(4,icur)<<"  ";
    
    cout<<endl<<" |     2000011 "<<setw(10)<<
      ( (mass(2000011) > 1e7) ? scientific : fixed)<<mass(2000011)<<fixed<<"  "; 
    for (int icur=1;icur<=8;icur++) cout<<setw(6)<<rvlmix(5,icur)<<"  ";
    
    cout<<endl<<" |     2000013 "<<setw(10)<<
      ( (mass(2000013) > 1e7) ? scientific : fixed)<<mass(2000013)<<fixed<<"  " ;
    for (int icur=1;icur<=8;icur++) cout<<setw(6)<<rvlmix(6,icur)<<"  ";
    
    cout<<endl<<" |     2000015 "<<setw(10)<<
      ( (mass(2000015) > 1e7) ? scientific : fixed)<<mass(2000015)<<fixed<<"  ";
    for (int icur=1;icur<=8;icur++) cout<<setw(6)<<rvlmix(7,icur)<<"  ";
  }
  cout<<endl;

  // Neutral scalars (sneutrinos)
  message(0,"","");

  // R-conserving:
  if (modsel(4) < 1) {
    cout<<" |  ~nu                 m";
    if (snumix.exists()) cout<<"   ~nu_e  ~nu_mu ~nu_tau";
    cout<<endl; 
    
    cout<<setprecision(3)<<" |     1000012 "<<setw(10)<<
      ( (mass(1000012) > 1e7) ? scientific : fixed)<<mass(1000012)<<fixed<<"  ";
    if (snumix.exists()) for (int icur=1;icur<=3;icur++) 
			   cout<<setw(6)<<snumix(1,icur)<<"  ";
    
    cout<<endl<<" |     1000014 "<<setw(10)<<
      ( (mass(1000014) > 1e7) ? scientific : fixed)<<mass(1000014)<<fixed<<"  ";
    if (snumix.exists()) for (int icur=1;icur<=3;icur++) 
			   cout<<setw(6)<<snumix(2,icur)<<"  ";
    
    cout<<endl<<" |     1000016 "<<setw(10)<<
      ( (mass(1000016) > 1e7) ? scientific : fixed)<<mass(1000016)<<fixed<<"  "; 
    if (snumix.exists()) for (int icur=1;icur<=3;icur++) 
			   cout<<setw(6)<<snumix(3,icur)<<"  ";
  }

  // R-violating
  else {
    cout<<" |  H0/~nu              m";
    if (snumix.exists()) cout<<"    H0_1    H0_2   ~nu_e  ~nu_mu ~nu_tau";
    cout<<endl; 
    
    cout<<setprecision(3)<<" |          25 "<<setw(10)<<
      ( (mass(25) > 1e7) ? scientific : fixed)<<mass(25)<<fixed<<"  ";
    if (rvhmix.exists()) for (int icur=1;icur<=5;icur++) 
			   cout<<setw(6)<<rvhmix(1,icur)<<"  ";
    
    cout<<endl<<" |          35 "<<setw(10)<<
      ( (mass(35) > 1e7) ? scientific : fixed)<<mass(35)<<fixed<<"  ";
    if (rvhmix.exists()) for (int icur=1;icur<=5;icur++) 
			   cout<<setw(6)<<rvhmix(2,icur)<<"  ";
    
    cout<<endl<<" |     1000012 "<<setw(10)<<
      ( (mass(1000012) > 1e7) ? scientific : fixed)<<mass(1000012)<<fixed<<"  ";
    if (rvhmix.exists()) for (int icur=1;icur<=5;icur++) 
			   cout<<setw(6)<<rvhmix(3,icur)<<"  ";
    
    cout<<endl<<" |     1000014 "<<setw(10)<<
      ( (mass(1000014) > 1e7) ? scientific : fixed)<<mass(1000014)<<fixed<<"  ";
    if (rvhmix.exists()) for (int icur=1;icur<=5;icur++) 
			   cout<<setw(6)<<rvhmix(4,icur)<<"  ";
    
    cout<<endl<<" |     1000016 "<<setw(10)<<
      ( (mass(1000016) > 1e7) ? scientific : fixed)<<mass(1000016)<<fixed<<"  "; 
    if (rvhmix.exists()) for (int icur=1;icur<=5;icur++) 
			   cout<<setw(6)<<rvhmix(5,icur)<<"  ";
  }
  cout<<endl;

  // Neutral pseudoscalars (RPV only)
  if (modsel(4) >= 1 && rvamix.exists()) {
    message(0,"","");
    cout<<" |  A0/~nu              m    A0_1    A0_2   ~nu_e  ~nu_mu ~nu_tau"<<endl; 

    cout<<setprecision(3)<<" |          36 "<<setw(10)<<
      ( (mass(36) > 1e7) ? scientific : fixed)<<mass(36)<<fixed<<"  ";
    for (int icur=1;icur<=5;icur++) cout<<setw(6)<<rvamix(1,icur)<<"  ";
    
    cout<<endl<<" |     1000017 "<<setw(10)<<
      ( (mass(1000017) > 1e7) ? scientific : fixed)<<mass(1000017)<<fixed<<"  ";
    for (int icur=1;icur<=5;icur++) cout<<setw(6)<<rvamix(2,icur)<<"  ";
    
    cout<<endl<<" |     1000018 "<<setw(10)<<
      ( (mass(1000018) > 1e7) ? scientific : fixed)<<mass(1000018)<<fixed<<"  ";
    for (int icur=1;icur<=5;icur++) cout<<setw(6)<<rvamix(3,icur)<<"  ";
    
    cout<<endl<<" |     1000019 "<<setw(10)<<
      ( (mass(1000019) > 1e7) ? scientific : fixed)<<mass(1000019)<<fixed<<"  ";
    for (int icur=1;icur<=5;icur++) cout<<setw(6)<<rvamix(4,icur)<<"  ";    
    cout<<endl;

  }

  // Neutral fermions (neutralinos)
  message(0,"","");

  // R-conserving:
  if (modsel(4) < 1) {
    cout<<" |  ~chi0               m      ~B    ~W_3    ~H_1    ~H_2"<<endl; 
    
    cout<<setprecision(3)<<" |     1000022 "<<setw(10)<<
      ( (mass(1000022) > 1e7) ? scientific : fixed)<<mass(1000022)<<fixed<<"  ";
    for (int icur=1;icur<=4;icur++) cout<<setw(6)<<nmix(1,icur)<<"  ";
    
    cout<<endl<<" |     1000023 "<<setw(10)<<
      ( (mass(1000023) > 1e7) ? scientific : fixed)<<mass(1000023)<<fixed<<"  ";
    for (int icur=1;icur<=4;icur++) cout<<setw(6)<<nmix(2,icur)<<"  ";
    
    cout<<endl<<" |     1000025 "<<setw(10)<<
      ( (mass(1000025) > 1e7) ? scientific : fixed)<<mass(1000025)<<fixed<<"  "; 
    for (int icur=1;icur<=4;icur++) cout<<setw(6)<<nmix(3,icur)<<"  ";
    
    cout<<endl<<" |     1000035 "<<setw(10)<<
      ( (mass(1000035) > 1e7) ? scientific : fixed)<<mass(1000035)<<fixed<<"  "; 
    for (int icur=1;icur<=4;icur++) cout<<setw(6)<<nmix(4,icur)<<"  ";
  }

  // R-violating
  else {
    cout<<" |  nu/~chi0            m    nu_e   nu_mu  nu_tau      ~B    ~W_3    ~H_1    ~H_2"<<endl; 
    
    cout<<setprecision(3)<<" |          12 "<<setw(10)<<
      ( (mass(12) > 1e7) ? scientific : fixed)<<mass(12)<<fixed<<"  ";
    for (int icur=1;icur<=7;icur++) cout<<setw(6)<<rvnmix(1,icur)<<"  ";
    
    cout<<endl<<" |          14 "<<setw(10)<<
      ( (mass(14) > 1e7) ? scientific : fixed)<<mass(14)<<fixed<<"  ";
    for (int icur=1;icur<=7;icur++) cout<<setw(6)<<rvnmix(2,icur)<<"  ";

    cout<<endl<<" |          16 "<<setw(10)<<
      ( (mass(16) > 1e7) ? scientific : fixed)<<mass(16)<<fixed<<"  ";
    for (int icur=1;icur<=7;icur++) cout<<setw(6)<<rvnmix(3,icur)<<"  ";

    cout<<endl<<" |     1000022 "<<setw(10)<<
      ( (mass(1000022) > 1e7) ? scientific : fixed)<<mass(1000022)<<fixed<<"  ";
    for (int icur=1;icur<=7;icur++) cout<<setw(6)<<rvnmix(4,icur)<<"  ";

    cout<<endl<<" |     1000023 "<<setw(10)<<
      ( (mass(1000023) > 1e7) ? scientific : fixed)<<mass(1000023)<<fixed<<"  ";
    for (int icur=1;icur<=7;icur++) cout<<setw(6)<<rvnmix(5,icur)<<"  ";
    
    cout<<endl<<" |     1000025 "<<setw(10)<<
      ( (mass(1000025) > 1e7) ? scientific : fixed)<<mass(1000025)<<fixed<<"  "; 
    for (int icur=1;icur<=7;icur++) cout<<setw(6)<<rvnmix(6,icur)<<"  ";
    
    cout<<endl<<" |     1000035 "<<setw(10)<<
      ( (mass(1000035) > 1e7) ? scientific : fixed)<<mass(1000035)<<fixed<<"  "; 
    for (int icur=1;icur<=7;icur++) cout<<setw(6)<<rvnmix(7,icur)<<"  ";
  }
  cout<<endl;

  // Charged fermions (charginos)
  message(0,"","");

  // R-conserving:
  if (modsel(4) < 1) {
    cout<<" |  ~chi+               m   U:   ~W      ~H  |  V:   ~W      ~H"
	<<endl; 
    
    cout<<setprecision(3)<<" |     1000024 "<<setw(10)<<
      ((mass(1000024) > 1e7) ? scientific : fixed)<<mass(1000024)<<fixed<<"    ";
    for (int icur=1;icur<=2;icur++) cout<<setw(6)<<umix(1,icur)<<"  ";
    cout<<"|   ";
    for (int icur=1;icur<=2;icur++) cout<<setw(6)<<vmix(1,icur)<<"  ";
    
    cout<<endl<<" |     1000037 "<<setw(10)<<
      ((mass(1000037) > 1e7) ? scientific : fixed)<<mass(1000037)<<fixed<<"    ";
    for (int icur=1;icur<=2;icur++) cout<<setw(6)<<umix(2,icur)<<"  ";
    cout<<"|   " ;
    for (int icur=1;icur<=2;icur++) cout<<setw(6)<<vmix(2,icur)<<"  ";
  }

  // R-violating
  else {
    cout<<" |  e+/~chi+            m   U:  eL+    muL+   tauL+     ~W+    ~H1+  |  V:  eR+    muR+   tauR+     ~W+    ~H2+"
	<<endl; 
    
    cout<<setprecision(3)<<" |         -11 "<<setw(10)<<
      ((mass(11) > 1e7) ? scientific : fixed)<<mass(11)<<fixed<<"    ";
    for (int icur=1;icur<=5;icur++) cout<<setw(6)<<rvumix(1,icur)<<"  ";
    cout<<"|   ";
    for (int icur=1;icur<=5;icur++) cout<<setw(6)<<rvvmix(1,icur)<<"  ";
    
    cout<<endl<<" |         -13 "<<setw(10)<<
      ((mass(13) > 1e7) ? scientific : fixed)<<mass(13)<<fixed<<"    ";
    for (int icur=1;icur<=5;icur++) cout<<setw(6)<<rvumix(2,icur)<<"  ";
    cout<<"|   " ;
    for (int icur=1;icur<=5;icur++) cout<<setw(6)<<rvvmix(2,icur)<<"  ";

    cout<<endl<<" |         -15 "<<setw(10)<<
      ((mass(15) > 1e7) ? scientific : fixed)<<mass(15)<<fixed<<"    ";
    for (int icur=1;icur<=5;icur++) cout<<setw(6)<<rvumix(3,icur)<<"  ";
    cout<<"|   " ;
    for (int icur=1;icur<=5;icur++) cout<<setw(6)<<rvvmix(3,icur)<<"  ";

    cout<<endl<<" |     1000024 "<<setw(10)<<
      ((mass(1000024) > 1e7) ? scientific : fixed)<<mass(1000024)<<fixed<<"    ";
    for (int icur=1;icur<=5;icur++) cout<<setw(6)<<rvumix(4,icur)<<"  ";
    cout<<"|   " ;
    for (int icur=1;icur<=5;icur++) cout<<setw(6)<<rvvmix(4,icur)<<"  ";

    cout<<endl<<" |     1000037 "<<setw(10)<<
      ((mass(1000037) > 1e7) ? scientific : fixed)<<mass(1000037)<<fixed<<"    ";
    for (int icur=1;icur<=5;icur++) cout<<setw(6)<<rvumix(5,icur)<<"  ";
    cout<<"|   " ;
    for (int icur=1;icur<=5;icur++) cout<<setw(6)<<rvvmix(5,icur)<<"  ";
  }
  cout<<endl;

  // Higgs bosons 
  message(0,"","");

  // R-conserving (R-violating case handled above, with sneutrinos)
  cout<<" |  Higgs      "<<endl;
  if (modsel(4) < 1) {
    cout<<setprecision(3);
    cout<<" |     alpha     ";
    if (alpha.exists()) cout<<setw(8)<<alpha()<<endl;
  }
  if (hmix.exists()) {
    cout<<" |     mu        ";
    if (hmix.exists(1)) cout<<setw(8)<<hmix(1)<<" (DRbar running value at Q = "<<hmix.q()<<" GeV)";
    else cout<<"  absent";
    cout<<"\n |     tan(beta) ";
    if (hmix.exists(2)) cout<<setw(8)<<hmix(2)<<" (DRbar running value at Q = "<<hmix.q()<<" GeV)";
    else cout<<"  absent";
    cout<<"\n |     v         ";
    if (hmix.exists(3)) cout<<setw(8)<<hmix(3)<<" (DRbar running value at Q = "<<hmix.q()<<" GeV)";
    else cout<<"  absent";
    cout<<"\n |     mA       ";
    if (hmix.exists(4)) cout<<setw(9)<<((abs(hmix(4)) > 1e5) ? scientific : fixed)<<hmix(4)
			    <<" (DRbar running value at Q = "<<fixed<<hmix.q()<<" GeV)";
    else cout<<"  absent";
    cout<<"\n";
  }

  // Gauge
  message(0,"","");
  if (gauge.exists()) {
    cout<<" |  Gauge      "<<endl;
    cout<<" |     g'        ";
    if (gauge.exists(1)) cout<<setw(8)<<gauge(1)<<" (DRbar running value at Q = "<<hmix.q()<<" GeV)";
    else cout<<"  absent";
    cout<<"\n |     g         ";
    if (gauge.exists(2)) cout<<setw(8)<<gauge(2)<<" (DRbar running value at Q = "<<hmix.q()<<" GeV)";
    else cout<<"  absent";
    cout<<"\n |     g3        ";
    if (gauge.exists(3)) cout<<setw(8)<<gauge(3)<<" (DRbar running value at Q = "<<hmix.q()<<" GeV)";
    else cout<<"  absent";
    cout<<"\n";
  }

  // Print footer  
  footerPrinted=false;
  message(0,"","");
  printFooter();
}

//--------------------------------------------------------------------------

// Check consistency of spectrum, unitarity of matrices, etc.

int SusyLesHouches::checkSpectrum() {

  if (! headerPrinted) printHeader();
  int ifail=0;
  bool foundModsel = modsel.exists();
  if (! foundModsel) {
    if (mass.exists()) return 1;
    else return 2;
  }

  // Step 1) Check MODSEL. Assign default values where applicable.
  if (!modsel.exists(1)) {
    message(1,"checkSpectrum","MODSEL(1) undefined. Assuming = 0",0);
    modsel.set(1,0);
    ifail=0;
  }
  if (!modsel.exists(3)) modsel.set(3,0);
  if (!modsel.exists(4)) modsel.set(4,0);
  if (!modsel.exists(5)) modsel.set(5,0);
  if (!modsel.exists(6)) modsel.set(6,0);
  if (!modsel.exists(11)) modsel.set(11,1);
  
  // Step 2) Check for existence / duplication of blocks

  //Global
  if (!minpar.exists()) {
      message(1,"checkSpectrum","MINPAR not found",0);
  }
  if (!sminputs.exists()) {
      message(1,"checkSpectrum","SMINPUTS not found",0);
  }
  if (!mass.exists()) {
      message(1,"checkSpectrum","MASS not found",0);
  }
  if (!gauge.exists()) {
      message(1,"checkSpectrum","GAUGE not found",0);
  }

  //SLHA1
  if (modsel(3) == 0 && modsel(4) == 0 && modsel(5) == 0 && modsel(6) == 0) {
    // Check for required SLHA1 blocks
    if (!staumix.exists() && !selmix.exists()) {
      message(1,"checkSpectrum","STAUMIX not found",0);
    };  
    if (!sbotmix.exists() && !dsqmix.exists()) {
      message(1,"checkSpectrum","SBOTMIX not found",0);
    };  
    if (!stopmix.exists() && !usqmix.exists()) {
      message(1,"checkSpectrum","STOPMIX not found",0);
    };  
    if (!nmix.exists()) {
      message(1,"checkSpectrum","NMIX not found",0);
    };  
    if (!umix.exists()) {
      message(1,"checkSpectrum","UMIX not found",0);
    };  
    if (!vmix.exists()) {
      message(1,"checkSpectrum","VMIX not found",0);
    };  
    if (!alpha.exists()) {
      message(1,"checkSpectrum","ALPHA not found",0);
    }
    if (!hmix.exists()) {
      message(1,"checkSpectrum","HMIX not found",0);
    }
    if (!msoft.exists()) {
      message(1,"checkSpectrum","MSOFT not found",0);
    }
  } 
  
  //RPV (+ FLV)
  else if (modsel(4) != 0) {
    // Check for required SLHA2 blocks (or see if can be extracted from SLHA1)
    if (!rvnmix.exists()) {
      if (nmix.exists()) {
	message(1,"checkSpectrum",
		"MODSEL 4 != 0 but NMIX given instead of RVNMIX",0);
	for (int i=1; i<=4; i++) {
	  if (i<=3) rvnmix.set(i,i,1.0);
	  for (int j=1; j<=4; j++)
	    rvnmix.set(i+3,j+3,nmix(i,j));	  
	} 
      } else {
	message(1,"checkSpectrum","MODSEL 4 != 0 but RVNMIX not found",0);
	ifail=-1;
      }
    }
    if (!rvumix.exists()) {
      if (umix.exists()) {
	message(1,"checkSpectrum",
		"MODSEL 4 != 0 but UMIX given instead of RVUMIX",0);
	for (int i=1; i<=3; i++) rvumix.set(i,i,1.0);	  
	for (int i=1; i<=2; i++) {
	  for (int j=1; j<=2; j++)
	    rvumix.set(i+3,j+3,umix(i,j));	  
	} 
      } else {
	message(1,"checkSpectrum","MODSEL 4 != 0 but RVUMIX not found",0);
	ifail=-1;
      }
    }
    if (!rvvmix.exists()) {
      if (vmix.exists()) {
	message(1,"checkSpectrum",
		"MODSEL 4 != 0 but VMIX given instead of RVVMIX",0);
	for (int i=1; i<=3; i++) rvvmix.set(i,i,1.0);	  
	for (int i=1; i<=2; i++) {
	  for (int j=1; j<=2; j++)
	    rvvmix.set(i+3,j+3,vmix(i,j));	  
	} 
      } else {
	message(1,"checkSpectrum","MODSEL 4 != 0 but RVVMIX not found",0);
	ifail=-1;
      }
    }
    if (!rvhmix.exists()) {
      if (alpha.exists()) {
	message(1,"checkSpectrum",
		"MODSEL 4 != 0 but ALPHA given instead of RVHMIX",0);
	rvhmix.set(1,1,cos(alpha()));
	rvhmix.set(1,2,sin(alpha()));
	rvhmix.set(2,1,-sin(alpha()));
	rvhmix.set(2,2,cos(alpha()));
	rvhmix.set(3,3,1.0);
	rvhmix.set(4,4,1.0);
	rvhmix.set(5,5,1.0);
      } else {
	message(1,"checkSpectrum","MODSEL 4 != 0 but RVHMIX not found",0);
	ifail=-1;
      }
    }
    if (!rvamix.exists()) {
      message(1,"checkSpectrum","MODSEL 4 != 0 but RVAMIX not found",0);
    }
    if (!rvlmix.exists()) {
      if (selmix.exists()) {
	message(1,"checkSpectrum",
		"MODSEL 4 != 0 but SELMIX given instead of RVLMIX",0);
	for (int i=1; i<=6; i++) {
	  for (int j=6; j<=6; j++)
	    rvlmix.set(i+1,j+2,selmix(i,j));	  
	} 	
      } if (staumix.exists()) {
	message(1,"checkSpectrum",
		"MODSEL 4 != 0 but STAUMIX given instead of RVLMIX",0);
	rvlmix.set(2,3,1.0);
	rvlmix.set(3,4,1.0);
	rvlmix.set(4,5,staumix(1,1));
	rvlmix.set(4,8,staumix(1,2));
	rvlmix.set(5,6,1.0);
	rvlmix.set(6,7,1.0);
	rvlmix.set(7,5,staumix(2,1));
	rvlmix.set(7,8,staumix(2,2));
      } else {
	message(1,"checkSpectrum","MODSEL 4 != 0 but RVLMIX not found",0);
	ifail=-1;
      }
    }
    if (!usqmix.exists()) {
      if (stopmix.exists()) {
	message(1,"checkSpectrum",
		"MODSEL 4 != 0 but STOPMIX given instead of USQMIX",0);
	usqmix.set(1,1, 1.0);
	usqmix.set(2,2, 1.0); 
	usqmix.set(4,4, 1.0);
	usqmix.set(5,5, 1.0);
	usqmix.set(3,3, stopmix(1,1));
	usqmix.set(3,6, stopmix(1,2));
	usqmix.set(6,3, stopmix(2,1));
	usqmix.set(6,6, stopmix(2,2));    
      } else {
	message(1,"checkSpectrum","MODSEL 4 != 0 but USQMIX not found",0);
	ifail=-1;
      }
    }
    if (!dsqmix.exists()) {
      if (sbotmix.exists()) {
	message(1,"checkSpectrum",
		"MODSEL 4 != 0 but SBOTMIX given instead of DSQMIX",0);
	dsqmix.set(1,1, 1.0);
	dsqmix.set(2,2, 1.0); 
	dsqmix.set(4,4, 1.0);
	dsqmix.set(5,5, 1.0);
	dsqmix.set(3,3, sbotmix(1,1));
	dsqmix.set(3,6, sbotmix(1,2));
	dsqmix.set(6,3, sbotmix(2,1));
	dsqmix.set(6,6, sbotmix(2,2));
      } else {
	message(1,"checkSpectrum","MODSEL 4 != 0 but DSQMIX not found",0);
	ifail=-1;
      }
    }
  }
  
  // FLV but not RPV (see above for FLV+RPV, below for FLV regardless of RPV)
  else if (modsel(6) != 0) {
    // Quark FLV
    if (modsel(6) != 2) {
      if (!usqmix.exists()) {
	message(1,"checkSpectrum","quark FLV on but USQMIX not found",0);
	ifail=-1;
      }
      if (!dsqmix.exists()) {
	message(1,"checkSpectrum","quark FLV on but DSQMIX not found",0);
	ifail=-1;
      }
    }
    // Lepton FLV
    if (modsel(6) != 1) {
      if (!upmns.exists()) {
	message(1,"checkSpectrum","lepton FLV on but UPMNSIN not found",0);
	ifail=-1;
      }
      if (!selmix.exists()) {
	message(1,"checkSpectrum","lepton FLV on but SELMIX not found",0);
	ifail=-1;
      }
      if (!snumix.exists() && !snsmix.exists()) {
	message(1,"checkSpectrum","lepton FLV on but SNUMIX not found",0);
	ifail=-1;
      }
    }
  }
  
  // CPV
  if (modsel(5) != 0) {
    if (!cvhmix.exists()) {
      message(1,"checkSpectrum","MODSEL 5 != 0 but CVHMIX not found",0);
      ifail=-1;
    }
  }
  
  // FLV (regardless of whether RPV or not)
  if (modsel(6) != 0) {
    // Quark FLV
    if (modsel(6) != 2) {
      if (!vckmin.exists()) {
	message(1,"checkSpectrum","quark FLV on but VCKMIN not found",0);
	ifail=-1;
      }
      if (!msq2in.exists()) {
	message(0,"checkSpectrum","note: quark FLV on but MSQ2IN not found",0);
	ifail=min(ifail,0);
      }
      if (!msu2in.exists()) {
	message(0,"checkSpectrum","note: quark FLV on but MSU2IN not found",0);
	ifail=min(ifail,0);
      }
      if (!msd2in.exists()) {
	message(0,"checkSpectrum","note: quark FLV on but MSD2IN not found",0);
	ifail=min(ifail,0);
      }
      if (!tuin.exists()) {
	message(0,"checkSpectrum","note: quark FLV on but TUIN not found",0);
	ifail=min(ifail,0);
      }
      if (!tdin.exists()) {
	message(0,"checkSpectrum","note: quark FLV on but TDIN not found",0);
	ifail=min(ifail,0);
      }
    }
    // Lepton FLV
    if (modsel(6) != 1) {
      if (!msl2in.exists()) {
	message(0,"checkSpectrum", 
                  "note: lepton FLV on but MSL2IN not found",0);
	ifail=min(ifail,0);
      }
      if (!mse2in.exists()) {
	message(0,"checkSpectrum",
                  "note: lepton FLV on but MSE2IN not found",0);
	ifail=min(ifail,0);
      }
      if (!tein.exists()) {
	message(0,"checkSpectrum",
                  "note: lepton FLV on but TEIN not found",0);
	ifail=min(ifail,0);
      }
    }
  }
  
  // Step 3) SLHA1 --> SLHA2 interoperability
  //Note: the mass basis is NOT mass-ordered in SLHA1, so be careful!
  //Here, the mass basis is hence by PDG code, not by mass-ordered value.

  if (stopmix.exists() && ! usqmix.exists() ) {
    //1000002 = ~uL, 1000004 = ~cL, 2000002 = ~uR, 2000004 = ~cR 
    usqmix.set(1,1, 1.0);
    usqmix.set(2,2, 1.0); 
    usqmix.set(4,4, 1.0);
    usqmix.set(5,5, 1.0);
    //Fill (1000006,2000006) sector from stopmix
    usqmix.set(3,3, stopmix(1,1));
    usqmix.set(3,6, stopmix(1,2));
    usqmix.set(6,3, stopmix(2,1));
    usqmix.set(6,6, stopmix(2,2));    
  };
  if (sbotmix.exists() && ! dsqmix.exists() ) {
    //1000001 = ~dL, 1000003 = ~sL, 2000001 = ~dR, 2000003 = ~sR 
    dsqmix.set(1,1, 1.0);
    dsqmix.set(2,2, 1.0); 
    dsqmix.set(4,4, 1.0);
    dsqmix.set(5,5, 1.0);
    //Fill (1000005,2000005) sector from sbotmix
    dsqmix.set(3,3, sbotmix(1,1));
    dsqmix.set(3,6, sbotmix(1,2));
    dsqmix.set(6,3, sbotmix(2,1));
    dsqmix.set(6,6, sbotmix(2,2));
  };
  if (staumix.exists() && ! selmix.exists() ) {
    //1000011 = ~eL, 1000013 = ~muL, 2000011 = ~eR, 2000013 = ~muR 
    selmix.set(1,1, 1.0);
    selmix.set(2,2, 1.0); 
    selmix.set(4,4, 1.0);
    selmix.set(5,5, 1.0);
    //Fill (1000015,2000015) sector from staumix
    selmix.set(3,3, staumix(1,1));
    selmix.set(3,6, staumix(1,2));
    selmix.set(6,3, staumix(2,1));
    selmix.set(6,6, staumix(2,2));
  };
  if (! snumix.exists() && ! snsmix.exists()) {
    //1000012 = ~nu_e, 1000014 = ~nu_mu, 1000016 = ~nu_tau
    snumix.set(1,1, 1.0);
    snumix.set(2,2, 1.0); 
    snumix.set(3,3, 1.0);
  };

  // Step 4) Check unitarity/orthogonality of mixing matrices

  //NMIX
  if (nmix.exists()) {
    for (int i=1;i<=4;i++) {
      double cn1=0.0;
      double cn2=0.0;
      for (int j=1;j<=4;j++) {
        cn1 += pow(nmix(i,j),2);
        cn2 += pow(nmix(j,i),2);
      }
      if (abs(1.0-cn1) > 1e-3 || abs(1.0-cn2) > 1e-3) { 
        ifail=2; 
        message(2,"checkSpectrum","NMIX is not unitary (wrong format?)",0);
	break;
      }
    }
  }

  //VMIX, UMIX
  if (vmix.exists() && umix.exists()) {
    // First check for non-standard "madgraph" convention
    // (2,2) entry not given explicitly
    for (int i=1;i<=2;i++) {
      double cu1=0.0;
      double cu2=0.0;
      double cv1=0.0;
      double cv2=0.0;
      for (int j=1;j<=2;j++) {
        cu1 += pow(umix(i,j),2);
        cu2 += pow(umix(j,i),2);
        cv1 += pow(vmix(i,j),2);
        cv2 += pow(vmix(j,i),2);
      }
      if (abs(1.0-cu1) > 1e-3 || abs(1.0-cu2) > 1e-3) { 
	cu1 += pow(umix(1,1),2);
	cu2 += pow(umix(1,1),2);
	if (abs(1.0-cu1) > 1e-3 || abs(1.0-cu2) > 1e-3) { 
	  ifail=max(1,ifail); 
	  message(2,"checkSpectrum","UMIX is not unitary (wrong format?)",0);
	  break;
	} else {
	  // Fix madgraph non-standard convention problem
	  message(1,"checkSpectrum","UMIX is not unitary (repaired)",0);
	  umix.set(2,2,umix(1,1));
	}
      }
      if (abs(1.0-cv1) > 1e-3 || abs(1.0-cv2) > 1e-3) { 
	cv1 += pow(vmix(1,1),2);
	cv2 += pow(vmix(1,1),2);
	if (abs(1.0-cv1) > 1e-3 || abs(1.0-cv2) > 1e-3) { 
	  ifail=max(1,ifail); 
	  message(2,"checkSpectrum","VMIX is not unitary (wrong format?)",0);
	  break;
	} else {
	  // Fix madgraph non-standard convention problem
	  message(1,"checkSpectrum","VMIX is not unitary (repaired)",0);
	  vmix.set(2,2,vmix(1,1));
	}
      }
    }
    
  }

  //STOPMIX, SBOTMIX
  if (stopmix.exists() && sbotmix.exists()) {
    for (int i=1;i<=2;i++) {
      double ct1=0.0;
      double ct2=0.0;
      double cb1=0.0;
      double cb2=0.0;
      for (int j=1;j<=2;j++) {
        ct1 += pow(stopmix(i,j),2);
        ct2 += pow(stopmix(j,i),2);
        cb1 += pow(sbotmix(i,j),2);
        cb2 += pow(sbotmix(j,i),2);
      }
      if (abs(1.0-ct1) > 1e-3 || abs(1.0-ct2) > 1e-3) { 
        ifail=-1; 
        message(2,"checkSpectrum","STOPMIX is not unitary (wrong format?)",0);
	break;
      }
      if (abs(1.0-cb1) > 1e-3 || abs(1.0-cb2) > 1e-3) { 
        ifail=-1; 
        message(2,"checkSpectrum","SBOTMIX is not unitary (wrong format?)",0);
	break;
      }
    }    
  }

  //STAUMIX
  if (staumix.exists()) {
    for (int i=1;i<=2;i++) {
      double ct1=0.0;
      double ct2=0.0;
      for (int j=1;j<=2;j++) {
        ct1 += pow(staumix(i,j),2);
        ct2 += pow(staumix(j,i),2);
      }
      if (abs(1.0-ct1) > 1e-3 || abs(1.0-ct2) > 1e-3) { 
        ifail=-1; 
        message(2,"checkSpectrum","STAUMIX is not unitary (wrong format?)",0);
	break;
      }
    }    
  }

  //DSQMIX
  if (dsqmix.exists()) {
    for (int i=1;i<=6;i++) {
      double sr=0.0;
      double sc=0.0;
      for (int j=1;j<=6;j++) {
        sr += pow(dsqmix(i,j),2);
        sc += pow(dsqmix(j,i),2);
      }
      if (abs(1.0-sr) > 1e-3 || abs(1.0-sc) > 1e-3) { 
        ifail=-1; 
        message(2,"checkSpectrum","DSQMIX is not unitary (wrong format?)",0);
	break;
      }
    }
  }

  //USQMIX
  if (usqmix.exists()) {
    for (int i=1;i<=6;i++) {
      double sr=0.0;
      double sc=0.0;
      for (int j=1;j<=6;j++) {
        sr += pow(usqmix(i,j),2);
        sc += pow(usqmix(j,i),2);
      }
      if (abs(1.0-sr) > 1e-3 || abs(1.0-sc) > 1e-3) { 
        ifail=-1; 
        message(2,"checkSpectrum","USQMIX is not unitary (wrong format?)",0);
	break;
      }
    }
  }

  //SELMIX
  if (selmix.exists()) {
    for (int i=1;i<=6;i++) {
      double sr=0.0;
      double sc=0.0;
      for (int j=1;j<=6;j++) {
        sr += pow(selmix(i,j),2);
        sc += pow(selmix(j,i),2);
      }
      if (abs(1.0-sr) > 1e-3 || abs(1.0-sc) > 1e-3) { 
        ifail=-1; 
        message(2,"checkSpectrum","SELMIX is not unitary (wrong format?)",0);
	break;
      }
    }
  }  //NMSSM:
  if (modsel(3) == 1) {
    //NMNMIX
    if ( nmnmix.exists() ) {
      for (int i=1;i<=5;i++) {
        double cn1=0.0;
        double cn2=0.0;
        for (int j=1;j<=5;j++) {
          cn1 += pow(nmnmix(i,j),2);
          cn2 += pow(nmnmix(j,i),2);
        }
        if (abs(1.0-cn1) > 1e-3 || abs(1.0-cn2) > 1e-3) { 
          ifail=-1;
          message(2,"checkSpectrum","NMNMIX is not unitary (wrong format?)",0);
	  break;
        }
      }
    }
    else {
      ifail=-1;
      message(1,"checkSpectrum","MODSEL 3 = 1 (NMSSM) but no NMNMIX found",0);
    }
    //NMAMIX
    if ( nmamix.exists() ) {
      for (int i=1;i<=2;i++) {
        double cn1=0.0;
        for (int j=1;j<=3;j++) {
          cn1 += pow(nmamix(i,j),2);
        }
        if (abs(1.0-cn1) > 1e-3) { 
          ifail=-1;
          message(2,"checkSpectrum","NMAMIX is not unitary (wrong format?)",0);
        }
      }
    }
    else {
      ifail=-1;
      message(1,"checkSpectrum","MODSEL 3 = 1 (NMSSM) but no NMAMIX found",0);
    }
    //NMHMIX
    if ( nmhmix.exists() ) {
      for (int i=1;i<=3;i++) {
        double cn1=0.0;
        double cn2=0.0;
        for (int j=1;j<=3;j++) {
          cn1 += pow(nmhmix(i,j),2);
          cn2 += pow(nmhmix(j,i),2);
        }
        if (abs(1.0-cn1) > 1e-3 || abs(1.0-cn2) > 1e-3) { 
          ifail=-1; 
          message(2,"checkSpectrum","NMHMIX is not unitary (wrong format?)",0);
        }
      }
    }
    else {
      ifail=-1;
      message(1,"checkSpectrum","MODSEL 3 = 1 (NMSSM) but no NMHMIX found",0);
    }
    //NMSSMRUN
    if (! nmssmrun.exists() ) {
      ifail=-1;
      message(2,"checkSpectrum","MODSEL 3 = 1 (NMSSM) but no NMSSMRUN found",
              0);
    }
  }
  
  //Check for documentation
  if (slhaRead && ! spinfo.exists(1)) spinfo.set(1,"unknown");
  if (slhaRead && ! spinfo.exists(2)) spinfo.set(2,"unknown");
  if (! slhaRead && ! spinfo.exists(1)) {
    spinfo.set(1,"DEFAULT");
    spinfo.set(2,"n/a");
  }

  //Give status 
  if (ifail >= 2) 
    message(0,"checkSpectrum","one or more serious problems were found");

  //Print Footer
  printFooter();

  //Return
  return ifail;
}

//--------------------------------------------------------------------------

// Check consistency of decay tables

int SusyLesHouches::checkDecays() {

  if (! headerPrinted) printHeader();
  int iFailDecays=0;
  
  // Loop over all particles read in
  for (int i = 0; i < int(decays.size()); ++i) { 
    
    // Shorthand
    LHdecayTable decTab = decays[i];
    int idRes = decTab.getId();
    double width = decTab.getWidth();
    if (width <= 0.0 || decTab.size() == 0) continue;
    
    // Check sum of branching ratios and phase spaces
    double sum = 0.0;
    double absSum = 0.0;
    int decSize = decTab.size();
    for (int j = 0; j < decSize; ++j) {
      
      double brat = decTab.getBrat(j);
      
      // Check phase space
      if (abs(brat) > 0.0) {
        vector<int> idDa = decTab.getIdDa(j);
        double massSum=abs(mass(idRes));      
        for (int k=0; k<int(idDa.size()); ++k) {
          if (mass.exists(idDa[k])) massSum -= mass(abs(idDa[k]));
          // If no MASS information read, use lowish values for check
          else if (abs(idDa[k]) == 24) massSum -=  79.0;
          else if (abs(idDa[k]) == 23) massSum -=  91.0;
          else if (abs(idDa[k]) ==  6) massSum -= 165.0;
          else if (abs(idDa[k]) ==  5) massSum -=   4.0;
          else if (abs(idDa[k]) ==  4) massSum -=   1.0;
        }
        if (massSum < 0.0) {
          // String containing decay name
          ostringstream errCode;
          errCode << idRes <<" ->";
          for (int jDa=0; jDa<int(idDa.size()); ++jDa) errCode<<" "<<idDa[jDa];
          message(1,"checkDecays",errCode.str()
                  +": Phase Space Closed, but BR != 0");
          iFailDecays = 1;
        }
        
      }

      // Sum up branching rations
      sum += brat;
      absSum += abs(brat);
      
    }

    if (abs(1.0-absSum) > 1e-6) {
      message(1,"checkDecays","sum(BR) != 1");
      cout << " | offending particle: "<<idRes<<" sum(BR) = "<<absSum<<endl;
      iFailDecays = 2;
    }
    
  }
  // End of loop over particles. Done.
  
  return iFailDecays;

}

//--------------------------------------------------------------------------

// Simple utility to print messages, warnings, and errors

void SusyLesHouches::message(int level, string place,string themessage,
  int line) {
  if (verbose == 0) return;
  //Send normal messages and warnings to stdout, errors to stderr.
  ostream* outstream = &cerr;
  if (level <= 1) outstream = &cout;
  // if (level == 2) { *outstream<<endl; }
  if (place != "") *outstream << " | (SLHA::"+place+") ";
  else *outstream << " | ";
  if (level == 1) *outstream<< "Warning: "; 
  if (level == 2) { *outstream <<"ERROR: "; } 
  if (line != 0) *outstream<< "line "<<line<<" - ";
  *outstream << themessage << endl;
  //  if (level == 2) *outstream <<endl;
  footerPrinted=false;
  return;
}

}

//==========================================================================





