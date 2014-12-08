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


## Comment.
class Comment:

  def __init__(self, lines, first_line, first_col, last_line, last_col, func):
    self.lines = lines
    self.first_line = first_line
    self.first_col = first_col
    self.last_line = last_line
    self.last_col = last_col
    self.func = func

  def __str__(self):
    return "<Comment for %s: [%d,%d:%d,%d] %s>" % (self.func, self.first_line, self.first_col, self.last_line, self.last_col, self.lines)


## Traverse the AST recursively starting from the current cursor.
#
#  @param cursor    A Clang parser cursor
#  @param comments  Array of comments found (of class Comment)
#  @param recursion Current recursion depth
def traverse_ast(cursor, comments, recursion=0):

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
    emit_comment = False

    comment = []
    comment_function = text
    comment_line_start = -1
    comment_line_end = -1
    comment_col_start = -1
    comment_col_end = -1

    for token in cursor.get_tokens():

      if token.cursor.kind == clang.cindex.CursorKind.COMPOUND_STMT:
        if not in_compound_stmt:
          in_compound_stmt = True
          expect_comment = True
          comment_line_end = -1
      else:
        if in_compound_stmt:
          in_compound_stmt = False
          emit_comment = True

      # tkind = str(token.kind)[str(token.kind).index('.')+1:]
      # ckind = str(token.cursor.kind)[str(token.cursor.kind).index('.')+1:]

      if in_compound_stmt:

        if expect_comment:

          extent = token.extent
          line_start = extent.start.line
          line_end = extent.end.line

          if token.kind == clang.cindex.TokenKind.PUNCTUATION and token.spelling == '{':
            pass

          elif token.kind == clang.cindex.TokenKind.COMMENT and (comment_line_end == -1 or (line_start == comment_line_end+1 and line_end-line_start == 0)):
            comment_line_end = line_end
            comment_col_end = extent.end.column

            if comment_line_start == -1:
              comment_line_start = line_start
              comment_col_start = extent.start.column
            comment.extend( token.spelling.split('\n') )

            # multiline comments are parsed in one go, therefore don't expect subsequent comments
            if line_end - line_start > 0:
              emit_comment = True
              expect_comment = False

          else:
            emit_comment = True
            expect_comment = False

      if emit_comment:

        comment = refactor_comment( comment )

        if len(comment) > 0:
          logging.debug("Comment found for function %s" % Colt(comment_function).magenta())
          comments.append( Comment(comment, comment_line_start, comment_col_start, comment_line_end, comment_col_end, comment_function) )

        comment = []
        comment_line_start = -1
        comment_line_end = -1
        comment_col_start = -1
        comment_col_end = -1

        emit_comment = False
        break

  else:

    logging.debug( "%s%s(%s)" % (indent, kind, text) )

  for child_cursor in cursor.get_children():
    traverse_ast(child_cursor, comments, recursion+1)


## Remove garbage from comments and convert special tags from THtml to Doxygen.
#
#  @param comment An array containing the lines of the original comment
def refactor_comment(comment):

  recomm = r'^(/{2,}|/\*)?\s*(.*?)\s*((/{2,})?\s*|\*/)$'

  new_comment = []
  insert_blank = False
  wait_first_non_blank = True
  for line_comment in comment:
    mcomm = re.search( recomm, line_comment )
    if mcomm:
      new_line_comment = mcomm.group(2)
      if new_line_comment == '':
        insert_blank = True
      else:
        if insert_blank and not wait_first_non_blank:
          new_comment.append('')
          insert_blank = False
        wait_first_non_blank = False
        new_comment.append( new_line_comment )
    else:
      assert False, 'Comment regexp does not match'

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

    logging.info('Input file: %s' % Colt(fn).magenta())
    index = clang.cindex.Index.create()
    translation_unit = index.parse(fn, args=['-x', 'c++'])

    comments = []
    traverse_ast( translation_unit.cursor, comments )
    for c in comments:
      logging.info("Comment found for %s:" % Colt(c.func).magenta())
      for l in c.lines:
        logging.info(
          Colt("[%d,%d:%d,%d] " % (c.first_line, c.first_col, c.last_line, c.last_col)).green() +
          "{%s}" % Colt(l).cyan()
        )

  return 0


if __name__ == '__main__':
  sys.exit( main( sys.argv[1:] ) )
