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
    print "%s%s(%s)" % (indent, Colt(kind).magenta(), Colt(text).blue())

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

          elif token.kind == clang.cindex.TokenKind.COMMENT and (last_comment_line == -1 or line_start == last_comment_line+1):
            #print Colt("%s  %s:%s = %s" % (indent, ckind, tkind, token.spelling)).green()
            last_comment_line = line_end
            print Colt("%s  [%d-%d]" % (indent, line_start, line_end)).cyan(),
            print Colt(token.spelling).yellow()

          else:
            expect_comment = False

      # else:
      #   print Colt("%s  %s:%s = %s" % (indent, ckind, tkind, token.spelling)).yellow()


  else:

    print "%s%s(%s)" % (indent, kind, text)

  # if cursor.kind == clang.cindex.CursorKind.CXX_METHOD:
  #   # here we are inside the method
  #   first_compound_stmt_found = False
  #   expect_comment_block = False
  #   print "%s ->->->->->->->->->->->->->->->->->->->->" % indent
  #   for token in cursor.get_tokens():

  #     print "%s %s::%s{%s}" % (indent, token.cursor.kind, token.kind, token.spelling)

  #     # find the first "compound statement"
  #     if token.cursor.kind == clang.cindex.CursorKind.COMPOUND_STMT:

  #       if not first_compound_stmt_found:
  #         # found first compound statement
  #         print ">>>> EXPECTING COMMENT ANY TIME NOW"
  #         first_compound_stmt_found = True
  #         expect_comment_block = True

  #     if expect_comment_block:
  #       if token.kind == clang.cindex.TokenKind.COMMENT:
  #         print "!!! [%s] %s" % (text, token.spelling)
  #       elif token.kind == clang.cindex.TokenKind.PUNCTUATION and token.spelling == '{':
  #         pass
  #       else:
  #         print "<<<< NOT EXPECTING COMMENT ANYMORE"
  #         expect_comment_block = False

  #   print "%s -<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<" % indent

  for child_cursor in cursor.get_children():
    traverse_ast(child_cursor, recursion+1)


## The main function.
#
#  **Note:** this program only has this function.
def main(argv):

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
    print Colt('[Error] Cannot find libclang, aborting').red()
    return 1

  # Loop over all files
  for fn in argv[1:]:

    index = clang.cindex.Index.create()
    translation_unit = index.parse(fn, args=['-x', 'c++'])
    traverse_ast( translation_unit.cursor )

  return 0


if __name__ == '__main__':
  sys.exit( main( sys.argv ) )
