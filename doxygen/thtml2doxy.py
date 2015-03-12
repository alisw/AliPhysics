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
#    `thtml2doxy [--stdout|-o] [-d] [--debug=DEBUG_LEVEL] file1 [file2 [file3...]]`
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
import hashlib
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
class Comment(object):

  def __init__(self, lines, first_line, first_col, last_line, last_col, indent, func, \
    append_empty=True):

    assert first_line > 0 and last_line >= first_line, 'Wrong line numbers'
    self.lines = lines
    self.first_line = first_line
    self.first_col = first_col
    self.last_line = last_line
    self.last_col = last_col
    self.indent = indent
    self.func = func
    self.append_empty = append_empty

  def has_comment(self, line):
    return line >= self.first_line and line <= self.last_line

  def __str__(self):
    return "<%s for %s: [%d,%d:%d,%d] %s>" % ( \
      self.__class__.__name__, self.func,
      self.first_line, self.first_col, self.last_line, self.last_col,
      self.lines)


## Prepend comment.
class PrependComment(Comment):

  def __init__(self, lines, first_line, first_col, last_line, last_col, indent, func, \
    append_empty=False):
    super(PrependComment, self).__init__( \
      lines, first_line, first_col, last_line, last_col, indent, func, append_empty)


## A data member comment.
class MemberComment:

  def __init__(self, text, comment_flag, array_size, first_line, first_col, func):
    assert first_line > 0, 'Wrong line number'
    assert comment_flag is None or comment_flag == '!' or comment_flag in [ '!', '||', '->' ]
    self.lines = [ text ]
    self.comment_flag = comment_flag
    self.array_size = array_size
    self.first_line = first_line
    self.first_col = first_col
    self.func = func

  def is_transient(self):
    return self.comment_flag == '!'

  def is_dontsplit(self):
    return self.comment_flag == '||'

  def is_ptr(self):
    return self.comment_flag == '->'

  def has_comment(self, line):
    return line == self.first_line

  def __str__(self):

    if self.is_transient():
      tt = '!transient! '
    elif self.is_dontsplit():
      tt = '!dontsplit! '
    elif self.is_ptr():
      tt = '!ptr! '
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

        comment = refactor_comment( comment, infilename=str(cursor.location.file) )

        if len(comment) > 0:
          logging.debug("Comment found for function %s" % Colt(comment_function).magenta())
          comments.append( Comment(comment, comment_line_start, comment_col_start, comment_line_end, comment_col_end, comment_indent, comment_function) )
        else:
          logging.debug('Empty comment found for function %s: collapsing' % Colt(comment_function).magenta())
          comments.append( Comment([''], comment_line_start, comment_col_start, comment_line_end, comment_col_end, comment_indent, comment_function) )
          #comments.append(RemoveComment(comment_line_start, comment_line_end))

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

  # Huge overkill: current line saved in "raw", previous in "prev"
  with open(str(cursor.location.file)) as fp:
    cur_line = 0
    for raw in fp:
      cur_line = cur_line + 1
      if cur_line == line_num:
        found = True
        break
      prev = raw

  assert found, 'A line that should exist was not found in file' % cursor.location.file

  recomm = r'(//(!|\|\||->)|///?)(\[(.+?)\])?<?!?\s*(.*?)\s*$'
  recomm_prevline = r'^\s*///\s*(.*?)\s*$'

  mcomm = re.search(recomm, raw)
  if mcomm:
    # If it does not match, we do not have a comment
    member_name = cursor.spelling;
    comment_flag = mcomm.group(2)  # !, ||, ->
    array_size = mcomm.group(4)
    text = mcomm.group(5)

    col_num = mcomm.start()+1;

    if array_size is not None and prev is not None:
      # ROOT arrays with comments already converted to Doxygen have the member description on the
      # previous line
      mcomm_prevline = re.search(recomm_prevline, prev)
      if mcomm_prevline:
        text = mcomm_prevline.group(1)
        comments.append(RemoveComment(line_num-1, line_num-1))

    logging.debug('Comment found for member %s' % Colt(member_name).magenta())

    comments.append( MemberComment(
      text,
      comment_flag,
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
#  Multi-line comments (/* ... */) that *immediately* follow a series of single-line comments
#  (*i.e.* without empty lines in-between) are also considered. A class description can eventually
#  be a series of single-line and multi-line comments, with no blank spaces between them, and always
#  starting with a single-line sequence.
#
#  The reason why they cannot start with a multi-line sequence is that those are commonly used to
#  display a copyright notice.
#
#  @param filename Name of the current file
#  @param comments Array of comments: new ones will be appended there
#  @param look_no_further_than_line Stop before reaching this line when looking for class comment
def comment_classdesc(filename, comments, look_no_further_than_line):

  # Single-line comment
  recomm = r'^\s*///?(\s*(.*?))\s*/*\s*$'

  # Multi-line comment (only either /* or */ on a single line)
  remlcomm_in = r'^\s*/\*\s*$'
  remlcomm_out = r'^\s*\*/\s*$'
  in_mlcomm = False

  reclass_doxy = r'(?i)^\s*\\(class|file):?\s*([^.]*)'
  class_name_doxy = None

  reauthor = r'(?i)^\s*\\?(authors?|origin):?\s*(.*?)\s*(,?\s*([0-9./-]+))?\s*$'
  redate = r'(?i)^\s*\\?date:?\s*([0-9./-]+)\s*$'
  rebrief = r'(?i)^\s*\\brief\s*(.*)\s*$'
  author = None
  date = None
  brief = None
  brief_len_threshold = 80

  comment_lines = []

  start_line = -1
  end_line = -1

  line_num = 0
  last_comm_line_num = 0

  is_macro = filename.endswith('.C')

  with open(filename, 'r') as fp:

    for raw in fp:

      line_num = line_num + 1

      if look_no_further_than_line is not None and line_num == look_no_further_than_line:
        logging.debug('Stopping at line %d while looking for class/file description' % \
          look_no_further_than_line)
        break

      if in_mlcomm == False and raw.strip() == '' and start_line > 0:
        # Skip empty lines
        continue

      stripped = strip_html(raw)
      mcomm = None
      this_comment = None

      if not in_mlcomm:
        mcomm = re.search(recomm, stripped)

      if last_comm_line_num == 0 or last_comm_line_num == line_num-1:

        if mcomm and not mcomm.group(2).startswith('#'):
          # Single-line comment
          this_comment = mcomm.group(1)
        elif start_line > -1:
          # Not a single-line comment. But it cannot be the first.
          if in_mlcomm == False:
            mmlcomm = re.search(remlcomm_in, stripped)
            if mmlcomm:
              in_mlcomm = True
              this_comment = ''
          else:
            mmlcomm = re.search(remlcomm_out, stripped)
            if mmlcomm:
              in_mlcomm = False
              this_comment = ''
            else:
              this_comment = stripped

      if this_comment is not None:

        if start_line == -1:

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

        end_line = line_num
        append = True

        mclass_doxy = re.search(reclass_doxy, this_comment)
        if mclass_doxy:
          class_name_doxy = mclass_doxy.group(2)
          append = False
        else:
          mauthor = re.search(reauthor, this_comment)
          if mauthor:
            author = mauthor.group(2)
            if date is None:
              # Date specified in the standalone \date field has priority
              date = mauthor.group(4)
            append = False
          else:
            mdate = re.search(redate, this_comment)
            if mdate:
              date = mdate.group(1)
              append = False
            else:
              mbrief = re.search(rebrief, this_comment)
              if mbrief:
                brief = mbrief.group(1)
                append = False

        if append:
          comment_lines.append( this_comment )

      else:
        if start_line > 0:
          break

      # This line had a valid comment
      last_comm_line_num = line_num

  if class_name_doxy is None:

    # No \class specified: guess it from file name
    reclass = r'^(.*/)?(.*?)(\..*)?$'
    mclass = re.search( reclass, filename )
    if mclass:
      class_name_doxy = mclass.group(2)
    else:
      assert False, 'Regexp unable to extract classname from file'

  # Macro or class?
  if is_macro:
    file_class_line = '\\file ' + class_name_doxy + '.C'
  else:
    file_class_line = '\\class ' + class_name_doxy

  if start_line > 0:

    prepend_to_comment = []

    # Prepend \class or \file specifier, then the \brief, then an empty line
    prepend_to_comment.append( file_class_line )

    if brief is not None:
      prepend_to_comment.append( '\\brief ' + brief )
    prepend_to_comment.append( '' )

    comment_lines = prepend_to_comment + comment_lines  # join lists

    # Append author and date if they exist
    if author is not None:
      comment_lines.append( '\\author ' + author )

    if date is not None:
      comment_lines.append( '\\date ' + date )

    # We should erase the "dumb" comments, such as "<class_name> class"
    comm_idx = 0
    regac = r'\s*%s\s+class\.?\s*' % class_name_doxy
    mgac = None
    for comm in comment_lines:
      mgac = re.search(regac, comm)
      if mgac:
        break
      comm_idx = comm_idx + 1
    if mgac:
      logging.debug('Removing dumb comment line: {%s}' % Colt(comment_lines[comm_idx]).magenta())
      del comment_lines[comm_idx]

    comment_lines = refactor_comment(comment_lines, do_strip_html=False, infilename=filename)

    # Now we look for a possible \brief
    if brief is None:
      comm_idx = 0
      for comm in comment_lines:
        if comm.startswith('\\class') or comm.startswith('\\file') or comm == '':
          pass
        else:
          if len(comm) <= brief_len_threshold:
            brief = comm
          break
        comm_idx = comm_idx + 1
      if brief is not None:
        comment_lines = refactor_comment(
          [ comment_lines[0], '\\brief ' + brief ] + comment_lines[1:comm_idx] + comment_lines[comm_idx+1:],
          do_strip_html=False, infilename=filename)

    logging.debug('Comment found for class %s' % Colt(class_name_doxy).magenta())
    comments.append(Comment(
      comment_lines,
      start_line, 1, end_line, 1,
      0, class_name_doxy
    ))

  else:

    logging.warning('No comment found for class %s: creating a dummy entry at the beginning' % \
      Colt(class_name_doxy).magenta())

    comments.append(PrependComment(
      [ file_class_line ],
      1, 1, 1, 1,
      0, class_name_doxy, append_empty=True
    ))


## Looks for a special ROOT ClassImp() entry.
#
#  Doxygen might get confused by `ClassImp()` entries as they are macros normally written without
#  the ending `;`: this function wraps the definition inside a condition in order to make Doxygen
#  ignore it.
#
#  @param filename Name of the current file
#  @param comments Array of comments: new ones will be appended there
def comment_classimp(filename, comments):

  recomm = r'^\s*///?(\s*.*?)\s*/*\s*$'

  line_num = 0
  reclassimp = r'^(\s*)Class(Imp|Def)\((.*?)\).*$'

  in_classimp_cond = False
  restartcond = r'^\s*///\s*\\cond\s+CLASSIMP\s*$'
  reendcond = r'^\s*///\s*\\endcond\s*$'

  with open(filename, 'r') as fp:

    # Array of tuples: classimp, startcond, endcond
    found_classimp = []

    # Reset to nothing found
    line_classimp = -1
    line_startcond = -1
    line_endcond = -1
    classimp_class = None
    classimp_indent = None

    for line in fp:

      line_num = line_num + 1

      mclassimp = re.search(reclassimp, line)
      if mclassimp:

        # Dump previous one if appropriate, and reset
        if line_classimp != -1:
          found_classimp.append( (line_classimp, line_startcond, line_endcond) )
          line_classimp = -1
          line_startcond = -1
          line_endcond = -1

        # Adjust indent
        classimp_indent = len( mclassimp.group(1) )

        line_classimp = line_num
        classimp_class = mclassimp.group(3)
        imp_or_def = mclassimp.group(2)
        logging.debug(
          'Comment found for ' +
          Colt( 'Class%s(' % imp_or_def ).magenta() +
          Colt( classimp_class ).cyan() +
          Colt( ')' ).magenta() )

      else:

        mstartcond = re.search(restartcond, line)
        if mstartcond:

          # Dump previous one if appropriate, and reset
          if line_classimp != -1:
            found_classimp.append( (line_classimp, line_startcond, line_endcond) )
            line_classimp = -1
            line_startcond = -1
            line_endcond = -1

          logging.debug('Found Doxygen opening condition for ClassImp')
          in_classimp_cond = True
          line_startcond = line_num

        elif in_classimp_cond:

          mendcond = re.search(reendcond, line)
          if mendcond:
            logging.debug('Found Doxygen closing condition for ClassImp')
            in_classimp_cond = False
            line_endcond = line_num

    # Dump previous one if appropriate, and reset (out of the loop)
    if line_classimp != -1:
      found_classimp.append( (line_classimp, line_startcond, line_endcond) )
      line_classimp = -1
      line_startcond = -1
      line_endcond = -1

    for line_classimp,line_startcond,line_endcond in found_classimp:

      # Loop over the ClassImp conditions we've found

      if line_startcond != -1:
        logging.debug('Looks like we are in a condition here %d,%d,%d' % (line_classimp, line_startcond, line_endcond))
        comments.append(Comment(
          ['\cond CLASSIMP'],
          line_startcond, 1, line_startcond, 1,
          classimp_indent, 'ClassImp/Def(%s)' % classimp_class,
          append_empty=False
        ))
      else:
        logging.debug('Looks like we are  NOT NOT  in a condition here %d,%d,%d' % (line_classimp, line_startcond, line_endcond))
        comments.append(PrependComment(
          ['\cond CLASSIMP'],
          line_classimp, 1, line_classimp, 1,
          classimp_indent, 'ClassImp/Def(%s)' % classimp_class
        ))

      if line_endcond != -1:
        comments.append(Comment(
          ['\endcond'],
          line_endcond, 1, line_endcond, 1,
          classimp_indent, 'ClassImp/Def(%s)' % classimp_class,
          append_empty=False
        ))
      else:
        comments.append(PrependComment(
          ['\endcond'],
          line_classimp+1, 1, line_classimp+1, 1,
          classimp_indent, 'ClassImp/Def(%s)' % classimp_class
        ))


## Traverse the AST recursively starting from the current cursor.
#
#  @param cursor    A Clang parser cursor
#  @param filename  Name of the current file
#  @param comments  Array of comments: new ones will be appended there
#  @param recursion Current recursion depth
#  @param in_func   True if we are inside a function or method
#  @param classdesc_line_limit  Do not look for comments after this line
#
#  @return A tuple containing the classdesc_line_limit as first item
def traverse_ast(cursor, filename, comments, recursion=0, in_func=False, classdesc_line_limit=None):

  # libclang traverses included files as well: we do not want this behavior
  if cursor.location.file is not None and str(cursor.location.file) != filename:
    logging.debug("Skipping processing of included %s" % cursor.location.file)
    return

  text = cursor.spelling or cursor.displayname
  kind = str(cursor.kind)[str(cursor.kind).index('.')+1:]

  is_macro = filename.endswith('.C')

  indent = ''
  for i in range(0, recursion):
    indent = indent + '  '

  if cursor.kind in [ clang.cindex.CursorKind.CXX_METHOD, clang.cindex.CursorKind.CONSTRUCTOR,
    clang.cindex.CursorKind.DESTRUCTOR, clang.cindex.CursorKind.FUNCTION_DECL ]:

    if classdesc_line_limit is None:
      classdesc_line_limit = cursor.location.line

    # cursor ran into a C++ method
    logging.debug( "%5d %s%s(%s)" % (cursor.location.line, indent, Colt(kind).magenta(), Colt(text).blue()) )
    comment_method(cursor, comments)
    in_func = True

  elif not is_macro and not in_func and \
    cursor.kind in [ clang.cindex.CursorKind.FIELD_DECL, clang.cindex.CursorKind.VAR_DECL ]:

    if classdesc_line_limit is None:
      classdesc_line_limit = cursor.location.line

    # cursor ran into a data member declaration
    logging.debug( "%5d %s%s(%s)" % (cursor.location.line, indent, Colt(kind).magenta(), Colt(text).blue()) )
    comment_datamember(cursor, comments)

  else:

    logging.debug( "%5d %s%s(%s)" % (cursor.location.line, indent, kind, text) )

  for child_cursor in cursor.get_children():
    classdesc_line_limit = traverse_ast(child_cursor, filename, comments, recursion+1, in_func, classdesc_line_limit)

  if recursion == 0:
    comment_classimp(filename, comments)
    comment_classdesc(filename, comments, classdesc_line_limit)

  return classdesc_line_limit


## Strip some HTML tags from the given string. Returns clean string.
#
#  @param s Input string
def strip_html(s):
  rehtml = r'(?i)</?(P|BR)/?>'
  return re.sub(rehtml, '', s)


## Remove garbage from comments and convert special tags from THtml to Doxygen.
#
#  @param comment An array containing the lines of the original comment
def refactor_comment(comment, do_strip_html=True, infilename=None):

  recomm = r'^(/{2,}|/\*)? ?(\s*)(.*?)\s*((/{2,})?\s*|\*/)$'
  regarbage = r'^(?i)\s*([\s*=_#-]+|(Begin|End)_Html)\s*$'

  # Support for LaTeX blocks spanning on multiple lines
  relatex = r'(?i)^((.*?)\s+)?(BEGIN|END)_LATEX([.,;:\s]+.*)?$'
  in_latex = False
  latex_block = False

  # Support for LaTeX blocks on a single line
  reinline_latex = r'(?i)(.*)BEGIN_LATEX\s+(.*?)\s+END_LATEX(.*)$'

  # Match <pre> (to turn it into the ~~~ Markdown syntax)
  reblock = r'(?i)^(\s*)</?PRE>\s*$'

  # Macro blocks for pictures generation
  in_macro = False
  current_macro = []
  remacro = r'(?i)^\s*(BEGIN|END)_MACRO(\((.*?)\))?\s*$'

  # Minimum indent level: scale back everything
  lowest_indent_level = None

  # Indentation threshold: if too much indented, don't indent at all
  indent_level_threshold = 7

  new_comment = []
  insert_blank = False
  wait_first_non_blank = True
  for line_comment in comment:

    # Strip some HTML tags
    if do_strip_html:
      line_comment = strip_html(line_comment)

    mcomm = re.search( recomm, line_comment )
    if mcomm:
      new_line_comment = mcomm.group(2) + mcomm.group(3)  # indent + comm

      # Check if we are in a macro block
      mmacro = re.search(remacro, new_line_comment)
      if mmacro:
        if in_macro:
          in_macro = False

          # Dump macro
          outimg = write_macro(infilename, current_macro) + '.png'
          current_macro = []

          # Insert image
          new_comment.append( '![Picture from ROOT macro](%s)' % (os.path.basename(outimg)) )

          logging.debug( 'Found macro for generating image %s' % Colt(outimg).magenta() )

        else:
          in_macro = True

        continue
      elif in_macro:
        current_macro.append( new_line_comment )
        continue

      mgarbage = re.search( regarbage, new_line_comment )

      if mgarbage is None and not mcomm.group(3).startswith('\\') and mcomm.group(3) != '':
        # not a special command line: count indent
        indent_level = len( mcomm.group(2) )
        if lowest_indent_level is None or indent_level < lowest_indent_level:
          lowest_indent_level = indent_level

        # if indentation level is too much, consider it zero
        if indent_level > indent_level_threshold:
          new_line_comment = mcomm.group(3)  # remove ALL indentation

      if new_line_comment == '' or mgarbage is not None:
        insert_blank = True
      else:
        if insert_blank and not wait_first_non_blank:
          new_comment.append('')
        insert_blank = False
        wait_first_non_blank = False

        # Postprocessing: LaTeX formulas in ROOT format
        # Marked by BEGIN_LATEX ... END_LATEX and they use # in place of \
        # There can be several ROOT LaTeX forumlas per line
        while True:
          minline_latex = re.search( reinline_latex, new_line_comment )
          if minline_latex:
            new_line_comment = '%s\\f$%s\\f$%s' % \
              ( minline_latex.group(1), minline_latex.group(2).replace('#', '\\'),
                minline_latex.group(3) )
          else:
            break

        # ROOT LaTeX: do we have a Begin/End_LaTeX block?
        # Note: the presence of LaTeX "closures" does not exclude the possibility to have a begin
        # block here left without a corresponding ending block
        mlatex = re.search( relatex, new_line_comment )
        if mlatex:

          # before and after parts have been already stripped
          l_before = mlatex.group(2)
          l_after = mlatex.group(4)
          is_begin = mlatex.group(3).upper() == 'BEGIN'  # if not, END

          if l_before is None:
            l_before = ''
          if l_after is None:
            l_after = ''

          if is_begin:

            # Begin of LaTeX part

            in_latex = True
            if l_before == '' and l_after == '':

              # Opening tag alone: mark the beginning of a block: \f[ ... \f]
              latex_block = True
              new_comment.append( '\\f[' )

            else:
              # Mark the beginning of inline: \f$ ... \f$
              latex_block = False
              new_comment.append(
                '%s \\f$%s' % ( l_before, l_after.replace('#', '\\') )
              )

          else:

            # End of LaTeX part
            in_latex = False

            if latex_block:

              # Closing a LaTeX block
              if l_before != '':
                new_comment.append( l_before.replace('#', '\\') )
              new_comment.append( '\\f]' )
              if l_after != '':
                new_comment.append( l_after )

            else:

              # Closing a LaTeX inline
              new_comment.append(
                '%s\\f$%s' % ( l_before.replace('#', '\\'), l_after )
              )

          # Prevent appending lines (we have already done that)
          new_line_comment = None

        # If we are not in a LaTeX block, look for <pre> tags and transform them into Doxygen code
        # blocks (using ~~~ ... ~~~). Only <pre> tags on a single line are supported
        if new_line_comment is not None and not in_latex:

          mblock = re.search( reblock, new_line_comment  )
          if mblock:
            new_comment.append( mblock.group(1)+'~~~' )
            new_line_comment = None

        if new_line_comment is not None:
          if in_latex:
            new_line_comment = new_line_comment.replace('#', '\\')
          new_comment.append( new_line_comment )

    else:
      assert False, 'Comment regexp does not match'

  # Fixing indentation level
  if lowest_indent_level is not None:
    logging.debug('Lowest indentation level found: %d' % lowest_indent_level)

    new_comment_indent = []
    reblankstart = r'^\s+'
    for line in new_comment:
      if re.search(reblankstart, line):
        new_comment_indent.append( line[lowest_indent_level:] )
      else:
        new_comment_indent.append( line )

    new_comment = new_comment_indent

  else:
    logging.debug('No indentation scaling applied')

  return new_comment


## Dumps an image-generating macro to the correct place. Returns a string with the image path,
#  without the extension.
#
#  @param infilename  File name of the source file
#  @param macro_lines Array of macro lines
def write_macro(infilename, macro_lines):

  # Calculate hash
  digh = hashlib.sha1()
  for l in macro_lines:
    digh.update(l)
    digh.update('\n')
  short_digest = digh.hexdigest()[0:7]

  infiledir = os.path.dirname(infilename)
  if infiledir == '':
    infiledir = '.'
  outdir = '%s/imgdoc' % infiledir
  outprefix = '%s/%s_%s' % (
    outdir,
    os.path.basename(infilename).replace('.', '_'),
    short_digest
  )
  outmacro = '%s.C' % outprefix

  # Make directory
  if not os.path.isdir(outdir):
    # do not catch: let everything die on error
    logging.debug('Creating directory %s' % Colt(outdir).magenta())
    os.mkdir(outdir)

  # Create file (do not catch errors either)
  with open(outmacro, 'w') as omfp:
    logging.debug('Writing macro %s' % Colt(outmacro).magenta())
    for l in macro_lines:
      omfp.write(l)
      omfp.write('\n')

  return outprefix


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
  restore_lines = None

  rindent = r'^(\s*)'

  def dump_comment_block(cmt, restore=None):
    text_indent = ''
    ask_skip_empty = False

    for i in range(0, cmt.indent):
      text_indent = text_indent + ' '

    for lc in cmt.lines:
      fhout.write('%s///' % text_indent )
      lc = lc.rstrip()
      if len(lc) != 0:
        fhout.write(' ')
        fhout.write(lc)
      fhout.write('\n')

    # Empty new line at the end of the comment
    if cmt.append_empty:
      fhout.write('\n')
      ask_skip_empty = True

    # Restore lines if possible
    if restore:
      for lr in restore:
        fhout.write(lr)
        fhout.write('\n')

    # Tell the caller whether it should skip the next empty line found
    return ask_skip_empty


  for line in fhin:

    line_num = line_num + 1

    # Find current comment
    prev_comm = comm
    comm = None
    comm_list = []
    for c in comments:
      if c.has_comment(line_num):
        comm = c
        comm_list.append(c)

    if len(comm_list) > 1:

      merged = True

      if len(comm_list) == 2:
        c1,c2 = comm_list
        if isinstance(c1, Comment) and isinstance(c2, Comment):
          c1.lines = c1.lines + c2.lines  # list merge
          comm = c1
          logging.debug('Two adjacent comments merged. Result: {%s}' % Colt(comm).cyan())
        else:
          merged = False
      else:
        merged = False

      if merged == False:
        logging.warning('Too many unmergeable comments on the same line (%d), picking the last one' % len(comm_list))
        for c in comm_list:
          logging.warning('>> %s' % c)
          comm = c  # considering the last one

    if comm:

      # First thing to check: are we in the same comment as before?
      if comm is not prev_comm and \
         isinstance(comm, Comment) and \
         isinstance(prev_comm, Comment) and \
         not isinstance(prev_comm, RemoveComment):

        # We are NOT in the same comment as before, and this comment is dumpable

        skip_empty = dump_comment_block(prev_comm, restore_lines)
        in_comment = False
        restore_lines = None
        prev_comm = None  # we have just dumped it: pretend it never existed in this loop

      #
      # Check type of comment and react accordingly
      #

      if isinstance(comm, MemberComment):

        # end comment block
        if in_comment:
          skip_empty = dump_comment_block(prev_comm, restore_lines)
          in_comment = False
          restore_lines = None

        non_comment = line[ 0:comm.first_col-1 ]

        if comm.array_size is not None or comm.is_dontsplit() or comm.is_ptr():

          # This is a special case: comment will be split in two lines: one before the comment for
          # Doxygen as "member description", and the other right after the comment on the same line
          # to be parsed by ROOT's C++ parser

          # Keep indent on the generated line of comment before member definition
          mindent = re.search(rindent, line)

          # Get correct comment flag, if any
          if comm.comment_flag is not None:
            cflag = comm.comment_flag
          else:
            cflag = ''

          # Get correct array size, if any
          if comm.array_size is not None:
            asize = '[%s]' % comm.array_size
          else:
            asize = ''

          # Write on two lines
          fhout.write('%s/// %s\n%s//%s%s\n' % (
            mindent.group(1),
            comm.lines[0],
            non_comment,
            cflag,
            asize
          ))

        else:

          # Single-line comments with the "transient" flag can be kept on one line in a way that
          # they are correctly interpreted by both ROOT and Doxygen

          if comm.is_transient():
            marker = '//!<!'  # compat with ROOT 5 and 6
          else:
            marker = '///<'

          fhout.write('%s%s %s\n' % (
            non_comment,
            marker,
            comm.lines[0]
          ))

      elif isinstance(comm, RemoveComment):
        # End comment block and skip this line
        if in_comment:
          skip_empty = dump_comment_block(prev_comm, restore_lines)
          in_comment = False
          restore_lines = None

      elif restore_lines is None:

        # Beginning of a new comment block of type Comment or PrependComment
        in_comment = True

        if isinstance(comm, PrependComment):
          # Prepare array of lines to dump right after the comment
          restore_lines = [ line.rstrip('\n') ]
          logging.debug('Commencing lines to restore: {%s}' % Colt(restore_lines[0]).cyan())
        else:
          # Extract the non-comment part and print it if it exists. If this is the first line of a
          # comment, it might happen something like `valid_code;  // this is a comment`.
          if comm.first_line == line_num:
            non_comment = line[ 0:comm.first_col-1 ].rstrip()
            if non_comment != '':
              fhout.write( non_comment + '\n' )

      elif isinstance(comm, Comment):

        if restore_lines is not None:
          # From the 2nd line on of comment to prepend
          restore_lines.append( line.rstrip('\n') )
          logging.debug('Appending lines to restore. All lines: {%s}' % Colt(restore_lines).cyan())

      else:
        assert False, 'Unhandled parser state: line=%d comm={%s} prev_comm={%s}' % \
          (line_num, comm, prev_comm)

    else:

      # Not a comment line

      if in_comment:

        # We have just exited a comment block of type Comment
        skip_empty = dump_comment_block(prev_comm, restore_lines)
        in_comment = False
        restore_lines = None

      # Dump the non-comment line
      line_out = line.rstrip('\n')
      if skip_empty:
        skip_empty = False
        if line_out.strip() != '':
          fhout.write( line_out + '\n' )
      else:
        fhout.write( line_out + '\n' )

  # Is there some comment left here?
  if restore_lines is not None:
    dump_comment_block(comm, restore_lines)

  # Is there some other comment beyond the last line?
  for c in comments:
    if c.has_comment(line_num+1):
      dump_comment_block(c, None)
      break


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
  include_flags = []
  try:
    opts, args = getopt.getopt( argv, 'odI:', [ 'debug=', 'stdout' ] )
    for o, a in opts:
      if o == '--debug':
        log_level = getattr( logging, a.upper(), None )
        if not isinstance(log_level, int):
          raise getopt.GetoptError('log level must be one of: DEBUG, INFO, WARNING, ERROR, CRITICAL')
      elif o == '-d':
        log_level = logging.DEBUG
      elif o == '-o' or o == '--stdout':
        output_on_stdout = True
      elif o == '-I':
        if os.path.isdir(a):
          include_flags.extend( [ '-I', a ] )
        else:
          logging.fatal('Include directory not found: %s' % Colt(a).magenta())
          return 2
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
    clang_args = [ '-x', 'c++' ]
    clang_args.extend( include_flags )
    translation_unit = index.parse(fn, args=clang_args)

    comments = []
    traverse_ast( translation_unit.cursor, fn, comments )
    for c in comments:

      logging.debug("Comment found for entity %s:" % Colt(c.func).magenta())

      if isinstance(c, MemberComment):

        if c.is_transient():
          flag_text = Colt('transient ').yellow()
        elif c.is_dontsplit():
          flag_text = Colt('dontsplit ').yellow()
        elif c.is_ptr():
          flag_text = Colt('ptr ').yellow()
        else:
          flag_text = ''

        if c.array_size is not None:
          array_text = Colt('arraysize=%s ' % c.array_size).yellow()
        else:
          array_text = ''

        logging.debug(
          "%s %s%s{%s}" % ( \
            Colt("[%d,%d]" % (c.first_line, c.first_col)).green(),
            flag_text,
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
