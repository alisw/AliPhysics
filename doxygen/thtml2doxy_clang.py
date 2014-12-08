#!/usr/bin/env python

## @package thtml2doxy_clang
#  Translates THtml C++ comments to Doxygen using libclang as parser.
#
#  This code relies on Python bindings for libclang: libclang's interface is pretty unstable, and
#  its Python bindings are unstable as well.
#
#  AST (Abstract Source Tree) traversal is performed entirely using libclang used as a C++ parser,
#  instead of attempting to write a parser ourselves.
#
#  This code (expecially AST traversal) was inspired by:
#
#   - [Implementing a code generator with libclang](http://szelei.me/code-generator/)
#     (this refers to API calls used here)
#   - [Parsing C++ in Python with Clang](http://eli.thegreenplace.net/2011/07/03/parsing-c-in-python-with-clang)
#     (outdated, API calls described there do not work anymore, but useful to understand some basic
#     concepts)
#
#  Usage:
#
#    `thtml2doxy_clang file1 [file2 [file3...]]`
#
#  @author Dario Berzano <dario.berzano@cern.ch>
#  @date 2014-12-05


import sys
import os
import re
import logging
import getopt
import clang.cindex


## Brain-dead color output for terminal.
class Colt(str):

  def red(self):
    return self.color('\033[31m')

  def green(self):
    return self.color('\033[32m')

  def yellow(self):
    return self.color('\033[33m')

  def blue(self):
    return self.color('\033[34m')

  def magenta(self):
    return self.color('\033[35m')

  def cyan(self):
    return self.color('\033[36m')

  def color(self, c):
    return c + self + '\033[m'


## Traverse the AST recursively starting from the current cursor.
#
#  @param cursor    A Clang parser cursor
#  @param recursion Current recursion depth
def traverse_ast(cursor, recursion=0):

  text = cursor.spelling or cursor.displayname
  kind = str(cursor.kind)[str(cursor.kind).index('.')+1:]

  indent = ''
  for i in range(0, recursion):
    indent = indent + '  '

  if cursor.kind == clang.cindex.CursorKind.CXX_METHOD:

    # cursor ran into a C++ method
    logging.debug( "%s%s(%s)" % (indent, Colt(kind).magenta(), Colt(text).blue()) )

    # we are looking for the following structure: method -> compound statement -> comment, i.e. we
    # need to extract the first comment in the compound statement composing the method

    in_compound_stmt = False
    expect_comment = False
    last_comment_line = -1

    for token in cursor.get_tokens():

      if token.cursor.kind == clang.cindex.CursorKind.COMPOUND_STMT:
        if not in_compound_stmt:
          in_compound_stmt = True
          expect_comment = True
          last_comment_line = -1
      else:
        if in_compound_stmt:
          in_compound_stmt = False
          break

      # tkind = str(token.kind)[str(token.kind).index('.')+1:]
      # ckind = str(token.cursor.kind)[str(token.cursor.kind).index('.')+1:]

      if in_compound_stmt:

        if expect_comment:

          extent = token.extent
          line_start = extent.start.line
          line_end = extent.end.line

          if token.kind == clang.cindex.TokenKind.PUNCTUATION and token.spelling == '{':
            pass

          elif token.kind == clang.cindex.TokenKind.COMMENT and (last_comment_line == -1 or (line_start == last_comment_line+1 and line_end-line_start == 0)):
            #print Colt("%s  %s:%s = %s" % (indent, ckind, tkind, token.spelling)).green()
            last_comment_line = line_end
            new_comment = refactor_comment(token.spelling)

            for comment_line in new_comment:
              logging.info(
                Colt("%s  [%d-%d]" % (indent, line_start, line_end)).green() +
                Colt(comment_line).cyan()
              )

            # multiline comments are parsed in one go, therefore don't expect subsequent comments
            if line_end - line_start > 0:
              expect_comment = False

          else:
            expect_comment = False

      # else:
      #   print Colt("%s  %s:%s = %s" % (indent, ckind, tkind, token.spelling)).yellow()


  else:

    logging.debug( "%s%s(%s)" % (indent, kind, text) )

  for child_cursor in cursor.get_children():
    traverse_ast(child_cursor, recursion+1)

## Remove garbage from comments and convert special tags from THtml to Doxygen.
#
#  @param comment The original comment
def refactor_comment(comment):

  resingle = r'^/{2,}\s*(.*?)\s*(/{2,})?\s*$'
  remulti_first = r'^/\*\s*(.*?)\s*\*?\s*$'
  remulti_last = r'^\s*(.*?)\s*\*/$'

  new_comment = comment.split('\n')

  if len(new_comment) == 1:
    msingle = re.search(resingle, comment)
    if msingle:
      new_comment[0] = msingle.group(1)

  else:

    for i in range(0, len(new_comment)):
      if i == 0:
        mmulti = re.search(remulti_first, new_comment[i])
        if mmulti:
          new_comment[i] = mmulti.group(1)
      elif i == len(new_comment)-1:
        mmulti = re.search(remulti_last, new_comment[i])
        if mmulti:
          new_comment[i] = mmulti.group(1)
      else:
        new_comment[i] = new_comment[i].strip()

  return new_comment


## The main function.
#
#  Return value is the executable's return value.
def main(argv):

  # Setup logging on stderr
  log_level = logging.WARNING
  logging.basicConfig(
    level=log_level,
    format='%(levelname)-8s %(funcName)-20s %(message)s',
    stream=sys.stderr
  )

  # Parse command-line options
  try:
    opts, args = getopt.getopt( argv, 'd', [ 'debug=' ] )
    for o, a in opts:
      if o == '--debug':
        log_level = getattr( logging, a.upper(), None )
        if not isinstance(log_level, int):
          raise getopt.GetoptError('log level must be one of: DEBUG, INFO, WARNING, ERROR, CRITICAL')
      elif o == '-d':
        log_level = logging.DEBUG
      else:
        assert False, 'Unhandled argument'
  except getopt.GetoptError as e:
    logging.fatal('Invalid arguments: %s' % e)
    return 1

  logging.getLogger('').setLevel(log_level)

  # Attempt to load libclang from a list of known locations
  libclang_locations = [
    '/usr/lib/llvm-3.5/lib/libclang.so.1',
    '/usr/lib/libclang.so',
    '/Library/Developer/CommandLineTools/usr/lib/libclang.dylib'
  ]
  libclang_found = False

  for lib in libclang_locations:
    if os.path.isfile(lib):
      clang.cindex.Config.set_library_file(lib)
      libclang_found = True
      break

  if not libclang_found:
    logging.fatal('Cannot find libclang')
    return 1

  # Loop over all files
  for fn in args:

    index = clang.cindex.Index.create()
    translation_unit = index.parse(fn, args=['-x', 'c++'])
    traverse_ast( translation_unit.cursor )

  return 0


if __name__ == '__main__':
  sys.exit( main( sys.argv[1:] ) )
