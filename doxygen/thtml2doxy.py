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
#    `thtml2doxy_clang [--stdout|-o] [-d] [--debug=DEBUG_LEVEL] file1 [file2 [file3...]]`
#
#  Parameters:
#
#   - `--stdout|-o`: output all on standard output instead of writing files in place
#   - `-d`: enable debug mode (very verbose output)
#   - `--debug=DEBUG_LEVEL`: set debug level to one of `DEBUG`, `INFO`, `WARNING`, `ERROR`,
#     `CRITICAL`
#
#  @author Dario Berzano, CERN
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

  def __init__(self, lines, first_line, first_col, last_line, last_col, indent, func):
    assert first_line > 0 and last_line >= first_line, 'Wrong line numbers'
    self.lines = lines
    self.first_line = first_line
    self.first_col = first_col
    self.last_line = last_line
    self.last_col = last_col
    self.indent = indent
    self.func = func

  def has_comment(self, line):
    return line >= self.first_line and line <= self.last_line

  def __str__(self):
    return "<Comment for %s: [%d,%d:%d,%d] %s>" % (self.func, self.first_line, self.first_col, self.last_line, self.last_col, self.lines)


## A data member comment.
class MemberComment:

  def __init__(self, text, is_transient, array_size, first_line, first_col, func):
    assert first_line > 0, 'Wrong line number'
    self.lines = [ text ]
    self.is_transient = is_transient
    self.array_size = array_size
    self.first_line = first_line
    self.first_col = first_col
    self.func = func

  def has_comment(self, line):
    return line == self.first_line

  def __str__(self):

    if self.is_transient:
      tt = '!transient! '
    else:
      tt = ''

    if self.array_size is not None:
      ars = '[%s] ' % self.array_size
    else:
      ars = ''

    return "<MemberComment for %s: [%d,%d] %s%s%s>" % (self.func, self.first_line, self.first_col, tt, ars, self.lines[0])


## A dummy comment that removes comment lines.
class RemoveComment(Comment):

  def __init__(self, first_line, last_line):
    assert first_line > 0 and last_line >= first_line, 'Wrong line numbers'
    self.first_line = first_line
    self.last_line = last_line
    self.func = '<remove>'

  def __str__(self):
    return "<RemoveComment: [%d,%d]>" % (self.first_line, self.last_line)


## Parses method comments.
#
#  @param cursor   Current libclang parser cursor
#  @param comments Array of comments: new ones will be appended there
def comment_method(cursor, comments):

  # we are looking for the following structure: method -> compound statement -> comment, i.e. we
  # need to extract the first comment in the compound statement composing the method

  in_compound_stmt = False
  expect_comment = False
  emit_comment = False

  comment = []
  comment_function = cursor.spelling or cursor.displayname
  comment_line_start = -1
  comment_line_end = -1
  comment_col_start = -1
  comment_col_end = -1
  comment_indent = -1

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

          if comment_indent == -1 or (extent.start.column-1) < comment_indent:
            comment_indent = extent.start.column-1

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

      if comment_line_start > 0:

        comment = refactor_comment( comment )

        if len(comment) > 0:
          logging.debug("Comment found for function %s" % Colt(comment_function).magenta())
          comments.append( Comment(comment, comment_line_start, comment_col_start, comment_line_end, comment_col_end, comment_indent, comment_function) )
        else:
          logging.debug('Empty comment for function %s marked for removal' % Colt(comment_function).magenta())
          comments.append(RemoveComment(comment_line_start, comment_line_end))

      else:
        logging.warning('No comment found for function %s' % Colt(comment_function).magenta())

      comment = []
      comment_line_start = -1
      comment_line_end = -1
      comment_col_start = -1
      comment_col_end = -1
      comment_indent = -1

      emit_comment = False
      break


## Parses comments to class data members.
#
#  @param cursor   Current libclang parser cursor
#  @param comments Array of comments: new ones will be appended there
def comment_datamember(cursor, comments):

  # Note: libclang 3.5 seems to have problems parsing a certain type of FIELD_DECL, so we revert
  # to a partial manual parsing. When parsing fails, the cursor's "extent" is not set properly,
  # returning a line range 0-0. We therefore make the not-so-absurd assumption that the datamember
  # definition is fully on one line, and we take the line number from cursor.location.

  line_num = cursor.location.line
  raw = None
  prev = None
  found = False

  # Huge overkill
  with open(str(cursor.location.file)) as fp:
    cur_line = 0
    for raw in fp:
      cur_line = cur_line + 1
      if cur_line == line_num:
        found = True
        break
      prev = raw

  assert found, 'A line that should exist was not found in file' % cursor.location.file

  recomm = r'(//(!)|///?)(\[(.*?)\])?<?\s*(.*?)\s*$'
  recomm_doxyary = r'^\s*///\s*(.*?)\s*$'

  mcomm = re.search(recomm, raw)
  if mcomm:
    # If it does not match, we do not have a comment
    member_name = cursor.spelling;
    is_transient = mcomm.group(2) is not None
    array_size = mcomm.group(4)
    text = mcomm.group(5)

    col_num = mcomm.start()+1;

    if array_size is not None and prev is not None:
      # ROOT arrays with comments already converted to Doxygen have the member description on the
      # previous line
      mcomm_doxyary = re.search(recomm_doxyary, prev)
      if mcomm_doxyary:
        text = mcomm_doxyary.group(1)
        comments.append(RemoveComment(line_num-1, line_num-1))

    logging.debug('Comment found for member %s' % Colt(member_name).magenta())

    comments.append( MemberComment(
      text,
      is_transient,
      array_size,
      line_num,
      col_num,
      member_name ))


## Parses class description (beginning of file).
#
#  The clang parser does not work in this case so we do it manually, but it is very simple: we keep
#  the first consecutive sequence of single-line comments (//) we find - provided that it occurs
#  before any other comment found so far in the file (the comments array is inspected to ensure
#  this).
#
#  Multi-line comments (/* ... */) are not considered as they are commonly used to display
#  copyright notice.
#
#  @param filename Name of the current file
#  @param comments Array of comments: new ones will be appended there
def comment_classdesc(filename, comments):

  recomm = r'^\s*///?(\s*.*?)\s*/*\s*$'

  reclass_doxy = r'(?i)^\s*\\class:?\s*(.*?)\s*$'
  class_name_doxy = None

  reauthor = r'(?i)^\s*\\?authors?:?\s*(.*?)\s*(,?\s*([0-9./-]+))?\s*$'
  redate = r'(?i)^\s*\\?date:?\s*([0-9./-]+)\s*$'
  author = None
  date = None

  comment_lines = []

  start_line = -1
  end_line = -1

  line_num = 0

  with open(filename, 'r') as fp:

    for raw in fp:

      line_num = line_num + 1

      if raw.strip() == '':
        # Skip empty lines
        end_line = line_num - 1
        continue

      stripped = strip_html(raw)
      mcomm = re.search(recomm, stripped)
      if mcomm:

        if start_line == -1 and len(comment_lines) == 0:

          # First line. Check that we do not overlap with other comments
          comment_overlaps = False
          for c in comments:
            if c.has_comment(line_num):
              comment_overlaps = True
              break

          if comment_overlaps:
            # No need to look for other comments
            break

          start_line = line_num

        append = True

        mclass_doxy = re.search(reclass_doxy, mcomm.group(1))
        if mclass_doxy:
          class_name_doxy = mclass_doxy.group(1)
          append = False
        else:
          mauthor = re.search(reauthor, mcomm.group(1))
          if mauthor:
            author = mauthor.group(1)
            if date is None:
              # Date specified in the standalone \date field has priority
              date = mauthor.group(2)
            append = False
          else:
            mdate = re.search(redate, mcomm.group(1))
            if mdate:
              date = mdate.group(1)
              append = False

        if append:
          comment_lines.append( mcomm.group(1) )

      else:
        if len(comment_lines) > 0:
          # End of our comment
          if end_line == -1:
            end_line = line_num - 1
          break

  if class_name_doxy is None:

    # No \class specified: guess it from file name
    reclass = r'^(.*/)?(.*?)(\..*)?$'
    mclass = re.search( reclass, filename )
    if mclass:
      class_name_doxy = mclass.group(2)
    else:
      assert False, 'Regexp unable to extract classname from file'

  # Prepend \class specifier (and an empty line)
  comment_lines[:0] = [ '\\class ' + class_name_doxy ]

  # Append author and date if they exist
  comment_lines.append('')

  if author is not None:
    comment_lines.append( '\\author ' + author )

  if date is not None:
    comment_lines.append( '\\date ' + date )

  comment_lines = refactor_comment(comment_lines, do_strip_html=False)
  logging.debug('Comment found for class %s' % Colt(class_name_doxy).magenta())
  comments.append(Comment(
    comment_lines,
    start_line, 1, end_line, 1,
    0, class_name_doxy
  ))


## Traverse the AST recursively starting from the current cursor.
#
#  @param cursor    A Clang parser cursor
#  @param filename  Name of the current file
#  @param comments  Array of comments: new ones will be appended there
#  @param recursion Current recursion depth
def traverse_ast(cursor, filename, comments, recursion=0):

  # libclang traverses included files as well: we do not want this behavior
  if cursor.location.file is not None and str(cursor.location.file) != filename:
    logging.debug("Skipping processing of included %s" % cursor.location.file)
    return

  text = cursor.spelling or cursor.displayname
  kind = str(cursor.kind)[str(cursor.kind).index('.')+1:]

  indent = ''
  for i in range(0, recursion):
    indent = indent + '  '

  if cursor.kind == clang.cindex.CursorKind.CXX_METHOD or cursor.kind == clang.cindex.CursorKind.CONSTRUCTOR or cursor.kind == clang.cindex.CursorKind.DESTRUCTOR:

    # cursor ran into a C++ method
    logging.debug( "%5d %s%s(%s)" % (cursor.location.line, indent, Colt(kind).magenta(), Colt(text).blue()) )
    comment_method(cursor, comments)

  elif cursor.kind == clang.cindex.CursorKind.FIELD_DECL:

    # cursor ran into a data member declaration
    logging.debug( "%5d %s%s(%s)" % (cursor.location.line, indent, Colt(kind).magenta(), Colt(text).blue()) )
    comment_datamember(cursor, comments)

  else:

    logging.debug( "%5d %s%s(%s)" % (cursor.location.line, indent, kind, text) )

  for child_cursor in cursor.get_children():
    traverse_ast(child_cursor, filename, comments, recursion+1)

  if recursion == 0:
    comment_classdesc(filename, comments)


## Strip some HTML tags from the given string. Returns clean string.
#
#  @param s Input string
def strip_html(s):
  rehtml = r'(?i)</?(P|H[0-9]|BR)/?>'
  return re.sub(rehtml, '', s)


## Remove garbage from comments and convert special tags from THtml to Doxygen.
#
#  @param comment An array containing the lines of the original comment
def refactor_comment(comment, do_strip_html=True):

  recomm = r'^(/{2,}|/\*)? ?(\s*.*?)\s*((/{2,})?\s*|\*/)$'
  regarbage = r'^(?i)\s*([\s*=-_#]+|(Begin|End)_Html)\s*$'

  new_comment = []
  insert_blank = False
  wait_first_non_blank = True
  for line_comment in comment:

    # Strip some HTML tags
    if do_strip_html:
      line_comment = strip_html(line_comment)

    mcomm = re.search( recomm, line_comment )
    if mcomm:
      new_line_comment = mcomm.group(2)
      mgarbage = re.search( regarbage, new_line_comment )

      if new_line_comment == '' or mgarbage is not None:
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


## Rewrites all comments from the given file handler.
#
#  @param fhin     The file handler to read from
#  @param fhout    The file handler to write to
#  @param comments Array of comments
def rewrite_comments(fhin, fhout, comments):

  line_num = 0
  in_comment = False
  skip_empty = False
  comm = None
  prev_comm = None

  rindent = r'^(\s*)'

  for line in fhin:

    line_num = line_num + 1

    # Find current comment
    prev_comm = comm
    comm = None
    for c in comments:
      if c.has_comment(line_num):
        comm = c

    if comm:

      if isinstance(comm, MemberComment):
        non_comment = line[ 0:comm.first_col-1 ]

        if comm.array_size is not None:

          mindent = re.search(rindent, line)
          if comm.is_transient:
            tt = '!'
          else:
            tt = ''

          # Special case: we need multiple lines not to confuse ROOT's C++ parser
          fhout.write('%s/// %s\n%s//%s[%s]\n' % (
            mindent.group(1),
            comm.lines[0],
            non_comment,
            tt,
            comm.array_size
          ))

        else:

          if comm.is_transient:
            tt = '!'
          else:
            tt = '/'

          fhout.write('%s//%s< %s\n' % (
            non_comment,
            tt,
            comm.lines[0]
          ))

      elif isinstance(comm, RemoveComment):
        # Do nothing: just skip line
        pass

      elif prev_comm is None:
        # Beginning of a new comment block of type Comment
        in_comment = True

        # Extract the non-comment part and print it if it exists
        non_comment = line[ 0:comm.first_col-1 ].rstrip()
        if non_comment != '':
          fhout.write( non_comment + '\n' )

    else:

      if in_comment:

        # We have just exited a comment block of type Comment
        in_comment = False

        # Dump revamped comment, if applicable
        text_indent = ''
        for i in range(0,prev_comm.indent):
          text_indent = text_indent + ' '

        for lc in prev_comm.lines:
          fhout.write( "%s/// %s\n" % (text_indent, lc) );
        fhout.write('\n')
        skip_empty = True

      line_out = line.rstrip('\n')
      if skip_empty:
        skip_empty = False
        if line_out.strip() != '':
          fhout.write( line_out + '\n' )
      else:
        fhout.write( line_out + '\n' )


## The main function.
#
#  Return value is the executable's return value.
def main(argv):

  # Setup logging on stderr
  log_level = logging.INFO
  logging.basicConfig(
    level=log_level,
    format='%(levelname)-8s %(funcName)-20s %(message)s',
    stream=sys.stderr
  )

  # Parse command-line options
  output_on_stdout = False
  try:
    opts, args = getopt.getopt( argv, 'od', [ 'debug=', 'stdout' ] )
    for o, a in opts:
      if o == '--debug':
        log_level = getattr( logging, a.upper(), None )
        if not isinstance(log_level, int):
          raise getopt.GetoptError('log level must be one of: DEBUG, INFO, WARNING, ERROR, CRITICAL')
      elif o == '-d':
        log_level = logging.DEBUG
      elif o == '-o' or o == '--stdout':
        logging.debug('Output on stdout instead of replacing original files')
        output_on_stdout = True
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
    traverse_ast( translation_unit.cursor, fn, comments )
    for c in comments:

      logging.debug("Comment found for entity %s:" % Colt(c.func).magenta())

      if isinstance(c, MemberComment):

        if c.is_transient:
          transient_text = Colt('transient ').yellow()
        else:
          transient_text = ''

        if c.array_size is not None:
          array_text = Colt('arraysize=%s ' % c.array_size).yellow()
        else:
          array_text = ''

        logging.debug(
          "%s %s%s{%s}" % ( \
            Colt("[%d,%d]" % (c.first_line, c.first_col)).green(),
            transient_text,
            array_text,
            Colt(c.lines[0]).cyan()
        ))

      elif isinstance(c, RemoveComment):

        logging.debug( Colt('[%d,%d]' % (c.first_line, c.last_line)).green() )

      else:
        for l in c.lines:
          logging.debug(
            Colt("[%d,%d:%d,%d] " % (c.first_line, c.first_col, c.last_line, c.last_col)).green() +
            "{%s}" % Colt(l).cyan()
          )

    try:

      if output_on_stdout:
        with open(fn, 'r') as fhin:
          rewrite_comments( fhin, sys.stdout, comments )
      else:
        fn_back = fn + '.thtml2doxy_backup'
        os.rename( fn, fn_back )

        with open(fn_back, 'r') as fhin, open(fn, 'w') as fhout:
          rewrite_comments( fhin, fhout, comments )

        os.remove( fn_back )
        logging.info("File %s converted to Doxygen: check differences before committing!" % Colt(fn).magenta())
    except (IOError,OSError) as e:
      logging.error('File operation failed: %s' % e)

  return 0


if __name__ == '__main__':
  sys.exit( main( sys.argv[1:] ) )
