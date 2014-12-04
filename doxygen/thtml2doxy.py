#!/usr/bin/env python

import sys
import os
import re

## @package thtml2doxy
#  Roughly translates THtml C++ comments to Doxygen.
#
#  Usage:
#
#    `thtml2doxy file1 [file2 [file3...]]`
#
#  @author Dario Berzano <dario.berzano@cern.ch>
#  @date 2014-12-04

## The main function.
#
#  **Note:** this program only has this function.
def main(argv):

  #reblock = r'^(\s*)//\s+(.*)\s*$'
  reblock = r'^//\s+(.*)\s*$'
  reignore = r'((Begin|End)_Html|^_+$)'
  restrip = r'(?i)</?P>|</?H[0-9]>|<BR/?>'
  refield = r'(?i)^(date|authors?):\s*(.*)$'
  reclass = r'/?([^./]+)\.[^.]+$'

  inblock = False
  endblock = False

  for fn in argv[1:]:

    mclass = re.search( reclass, fn )
    classname = mclass.group(1)

    fnout = fn + '.doxyout'
    try:
      with open(fn, 'r') as fp, open(fnout, 'w') as fout:
        print '%s...' % fn
        for rawline in fp:

          mblock = re.match(reblock, rawline)

          if mblock and not endblock and not inblock:
            # beginning of a block
            fout.write('/// \class %s\n' % classname)
            inblock = True
          elif not mblock and inblock:
            # end of a block
            inblock = False
            endblock = True

          if inblock:

            #indent = mblock.group(1)
            #text = mblock.group(2)
            indent = ''
            text = mblock.group(1)

            if re.match(reignore, text):
              pass
            else:

              stripped = re.sub( restrip, '', text ).strip()

              mfield = re.match(refield, stripped)
              if mfield:
                field = mfield.group(1).lower()
                value = mfield.group(2)
                fout.write( '%s/// \%s: %s\n' % (indent, field, value) )

              else:
                # print by restoring indent and with correct comment
                if len(stripped) > 0:
                  stripped = ' ' + stripped
                fout.write( '%s///%s\n' % (indent, stripped) )

          else:
            # no comment
            fout.write( rawline )

    except IOError as e:
      print '[ERROR] cannot open %s: %s' % (fn, e)

    try:
      os.remove( fn )
      os.rename( fnout, fn )
    except OSError as e:
      print '[ERROR] renaming/removing: %s (aborting)' % e
      return 1

  return 0

if __name__ == '__main__':
  sys.exit( main( sys.argv ) )
