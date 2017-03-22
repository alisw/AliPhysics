///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2011
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
// $Rev:: 208                       $: revision of last commit
// $Author:: jnystrand          $: author of last commit
// $Date:: 2015-08-09 22:29:36 +#$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef INPUTPARSER_H
#define INPUTPARSER_H

#include <string>
#include <typeinfo>
#include <iostream>
#include <map>

#include <reportingUtils.h>

class inputParser
{
public:
  
  /** Constructor */
  inputParser();

  /** Destructor */
  ~inputParser();

  /** Parse a file */
  int parseFile(std::string filename);

  /** Parse a file */
  int parseString(std::string str);
  
  /** Add parameter to pass */
  void addIntParameter(std::string name, int *var, bool required = true);

  /** Add parameter to pass */
  void addUintParameter(std::string name, unsigned int *var, bool required = true);

  /** Add parameter to pass */
  void addFloatParameter(std::string name, float *var, bool required = true);

  /** Add parameter to pass */
  void addDoubleParameter(std::string name, double *var, bool required = true);

  /** Add parameter to pass */
  void addBoolParameter(std::string name, bool *var, bool required = true);
 
  /** Add parameter to pass */
  void addStringParameter(std::string name, std::string *var, bool required = true);
  
  /** Print info */
  void printParameterInfo(std::ostream &out = std::cout);
  
  /** Validate */
  bool validateParameters(std::ostream &errOut = std::cerr);
  
  /** Add a parameter */
  template<typename S>
  inline void addParameter(S &param);
  
  /** Add a parameter */
  template<typename P>
  inline void addParameter(const std::string &name, P *varPtr, bool required = false);

private:
  
  template <class T>
  class _parameter
  {
  public:
    _parameter(std::string name, T *val, bool required = true, bool found = false) : _name(name), _val(val), _required(required), _found(found){}
    
    bool operator==(const _parameter &rhs) const { return _name == rhs._name; }
    
    bool operator<(const _parameter &rhs) const { return _name.c_str()[0] < rhs._name.c_str()[0]; }
    
    void printParameterInfo(std::ostream &out = std::cout) 
    {
      out << std::boolalpha << _name << "\t\t";
      if(_found)
      {
	out << *_val << std::endl;
      }
      else
      {
	out << "NOT FOUND" << std::endl;
      }
      out << std::noboolalpha;
    }
    
    
    std::string _name;
    T *_val;
    bool _required;
    bool _found;
  };
  
  std::map<std::string, _parameter<int> > _intParameters;
  std::map<std::string, _parameter<unsigned int> > _uintParameters;
  std::map<std::string, _parameter<float> > _floatParameters;
  std::map<std::string, _parameter<double> > _doubleParameters;
  std::map<std::string, _parameter<bool> > _boolParameters;
  std::map<std::string, _parameter<std::string> > _stringParameters;
  
};

template<typename S>
void inputParser::addParameter(S& param)
{
  addParameter(param.name(), param.ptr(), param.required());

}

template<typename P>
void inputParser::addParameter(const std::string& name, P* /*varPtr*/, bool /*required*/)
{
  printWarn << "Trying to add unknown parameter type with name: " << name;
}


template<>
inline void inputParser::addParameter(const std::string& name, int * varPtr, bool required)
{
  addIntParameter(name, varPtr, required);
}

template<>
inline void inputParser::addParameter(const std::string& name, unsigned int * varPtr, bool required)
{
  addUintParameter(name, varPtr, required);
}

template<>
inline void inputParser::addParameter(const std::string& name, float * varPtr, bool required)
{
  addFloatParameter(name, varPtr, required);
}

template<>
inline void inputParser::addParameter(const std::string& name, double * varPtr, bool required)
{
  addDoubleParameter(name, varPtr, required);
}

template<>
inline void inputParser::addParameter(const std::string& name, bool * varPtr, bool required)
{
  addBoolParameter(name, varPtr, required);
}

template<>
inline void inputParser::addParameter(const std::string& name, std::string * varPtr, bool required)
{
  addStringParameter(name, varPtr, required);
}

#endif // INPUTPARSER_H
