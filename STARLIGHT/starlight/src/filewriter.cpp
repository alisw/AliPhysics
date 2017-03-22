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


#include <iostream>
#include <exception>
#include <cstdlib>

#include "filewriter.h"


using namespace std;


fileWriter::fileWriter()
: _fileName(""),
  _fileStream()
{ }


fileWriter::fileWriter(const string& fileName) :
        _fileName(fileName)
        ,_fileStream(fileName.c_str())
{ }


fileWriter::~fileWriter()
{ }


int fileWriter::open()
{
    try
    {
        _fileStream.open(_fileName.c_str());
    }
    catch (const ios::failure & error)
    {
        cerr << "I/O exception: " << error.what() << endl;
        return EXIT_FAILURE;
    }
    return 0;
}


int fileWriter::open(const string& fileName)
{
    _fileName = fileName;
    return open();
}


int fileWriter::close()
{
    try
    {
        _fileStream.close();
    }
    catch (const ios::failure & error)
    {
        cerr << "I/O exception: " << error.what() << endl;
        return EXIT_FAILURE;
    }
    return 0;
}
