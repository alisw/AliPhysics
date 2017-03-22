///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of starlight.
//
//    starlight is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    starlight is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with starlight. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//
// File and Version Information:
// $Rev:: 28                          $: revision of last commit
// $Author:: bgrube                   $: author of last commit
// $Date:: 2010-12-10 19:30:01 +0100 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef FILEWRITER_H
#define FILEWRITER_H


#include <string>
#include <fstream>

#include "upcevent.h"


class fileWriter
{
public:

    /** Default constructor */
    fileWriter();

    /** Constructor with filename */
    fileWriter(const std::string& fileName);

    /** Destructor */
    virtual ~fileWriter();

    /** open the file */
    int open();

    /** open file with given filename */
    int open(const std::string& fileName);
    
    /** close the file */
    int close();

    /** Set the filename we're writing to */
    void setFileName(const std::string& fileName) { _fileName = fileName; }

protected:

   /** The file name */
    std::string _fileName;

    /** The file stream */ 
    std::ofstream _fileStream;
};


#endif  // FILEWRITER_H
