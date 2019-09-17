<html>
<head>
<title>The Settings Scheme</title>
<link rel="stylesheet" type="text/css" href="pythia.css"/>
<link rel="shortcut icon" href="pythia32.gif"/>
</head>
<body>

<script language=javascript type=text/javascript>
function stopRKey(evt) {
var evt = (evt) ? evt : ((event) ? event : null);
var node = (evt.target) ? evt.target :((evt.srcElement) ? evt.srcElement : null);
if ((evt.keyCode == 13) && (node.type=="text"))
{return false;}
}

document.onkeypress = stopRKey;
</script>
<?php
if($_POST['saved'] == 1) {
if($_POST['filepath'] != "files/") {
echo "<font color='red'>SETTINGS SAVED TO FILE</font><br/><br/>"; }
else {
echo "<font color='red'>NO FILE SELECTED YET.. PLEASE DO SO </font><a href='SaveSettings.php'>HERE</a><br/><br/>"; }
}
?>

<form method='post' action='SettingsScheme.php'>
<h2>The Settings Scheme</h2> 
<ol id="toc">
  <li><a href="#section0">Concepts</a></li>
  <li><a href="#section1">Operation</a></li>
  <li><a href="#section2">The public methods</a></li>
</ol>

 
The <code>Settings</code> class keeps track of all the flags, modes, 
parameters and words used during the event generation. As such, it 
serves all the <code>Pythia</code> program elements from one central 
repository. Accessing it allows the user to modify the generator 
behaviour. 
 
<p/> 
Each <code>Pythia</code> object has a public member <code>settings</code> 
of the <code>Settings</code> class. Therefore you access the 
settings methods as <code>pythia.settings.command(argument)</code>, 
assuming that <code>pythia</code> is an instance of the <code>Pythia</code> 
class. Further, for the most frequent user tasks, <code>Pythia</code> 
methods have been defined, so that <code>pythia.command(argument)</code> 
would work, see further below. 
 
<p/> 
The central section on this page is the Operation one. The preceding 
concepts section is there mainly to introduce the basic structure and 
the set of properties that can be accessed. The subsequent sections 
provide a complete listing of the existing public methods, which most 
users probably will have little interaction with. 
 
<a name="section0"></a> 
<h3>Concepts</h3> 
 
We distinguish eight kinds of user-modifiable variables, by the way 
they have to be stored: 
<ol> 
<li>Flags are on/off switches, and are stored as <code>bool</code>.</li> 
<li>Modes corresponds to a finite enumeration of separate options, 
and are stored as <code>int</code>.</li> 
<li>Parameters take a continuum of values, and are stored as 
<code>double</code>. The shorthand notation parm is used in the C++ 
code and XML tags.</li> 
<li>Words are simple character strings and are stored as 
<code>string</code>. No double quotation marks &quot; or braces { } 
may appear inside a word, and commas , will take a spcial role next 
so should also be avoided. Normally the input string is expected not to 
contain any blanks or equal signs, but if it does it must be enclosed 
in braces { }.</li> 
<li>Vectors of flags take a variable length, and are stored as 
<code>vector&lt;bool&gt;</code>. The shorthand notation fvec is used 
in the C++ code and XML tags. When the vector is input as a string 
it should be given as a comma-separated list, either containing 
no blanks or else enclosed in braces { }.</li> 
<li>Vectors of modes take a variable length, and are stored as 
<code>vector&lt;int&gt;</code>. The shorthand notation mvec is used 
in the C++ code and XML tags. When the vector is input as a string 
it should be given as a comma-separated list, either containing 
no blanks or else enclosed in braces { }.</li> 
<li>Vectors of parameters take a variable length and for each element 
a continuum of values, and are stored as <code>vector&lt;double&gt;</code>. 
The shorthand notation pvec is used in the C++ code and XML tag. 
When the vector is input as a string it should be given as a comma-separated 
list, either containing no blanks or else enclosed in braces { }.</li> 
<li>Vectors of words take a variable length, and are stored as 
<code>vector&lt;string&gt;</code>. The shorthand notation wvec is used 
in the C++ code and XML tags. When the vector is input as a string 
it should be given as a comma-separated list, either containing 
no blanks or else enclosed in braces { }.</li> 
</ol> 
Note the special role played by the braces { } to enclose words or lists 
that are allowed to contain blanks and equal signs, and of commas , to 
separate the fields of the list, in analogy with how C++ arrays can be 
initialized. You should not be using these three characters for any other 
purposes. Input of a vector can be split across several lines, until a 
close brace } is found that matches the open brace {. If no such is found 
the program will abort, so beware. The double quotation mark &quot; is 
avoided since it is already used for other purposes. Also note that all 
shorthands have been chosen four letters long. 
<p/> 
In general, each variable stored in <code>Settings</code> is associated 
with four kinds of information: 
<ul> 
<li>The variable name, of the form <code>class:name</code> 
(or <code>file:name</code>, usually these agree), e.g. 
<code>TimeShower:pTmin</code>. The class/file part usually identifies 
the <code>.xml</code> file where the variable is defined, and the part of 
the program where it is used, but such a connection cannot be strictly 
upheld, since e.g. the same variable may be used in a few different 
cases (even if most of them are not).</li> 
<li>The default value, set in the original declaration, and intended 
to represent a reasonable choice.</li> 
<li>The current value, which differs from the default when the user so 
requests.</li> 
<li>An allowed range of values, represented by meaningful 
minimum and maximum values. This has no sense for a <code>flag</code>, 
an <code>fvec</code> a <code>word</code> or a <code>wvec</code> 
(and is not used there), is 
usually rather well-defined for a <code>mode</code> or <code>mvec</code>, 
but less so for a <code>parm</code> or <code>pvec</code>. Often the allowed 
range exaggerates the degree of our current knowledge, so as not to restrict 
too much what the user can do. One may choose not to set the lower or 
upper limit, in which case the range is open-ended. 
<br/>Normally input values outside the allowed range are changed 
to the be at nearest limit. For <code>mode</code>s only, a further 
boolean is stored to tell whether this should be allowed, or whether 
out-of-range inputs should be forbidden, to the extent that the whole 
PYTHIA initialization would abort. The latter applies to those modes 
that have been defined with the <code>modepick</code> label in the 
<code>xmldoc/*.xml</code> files, and where maximal and minimal values 
have been specified. Such labels are used to represent a discrete set 
of options, and so any value outside the allowed range is just plain 
wrong. Also attempts to change <code>modefix</code> fixed-value 
<code>mode</code>s lead to aborts. By contrast those defined with 
<code>mode</code> or <code>modeopen</code> follow the normal rules 
of being reset to fall into the allowed range, without any warnings.</li> 
</ul> 
 
<p/> 
Technically, the <code>Settings</code> class is implemented with the 
help of eight separate maps, one for each kind of variable, with the 
variable <code>name</code> used as key. 
 
<a name="section1"></a> 
<h3>Operation</h3> 
 
The normal flow of setting values is: 
 
<ol> 
 
<p/> <li> 
When a <code>Pythia</code> object <code>pythia </code>is created, 
the member <code>pythia.settings</code> is asked to scan the files 
listed in the <code>Index.xml</code> file in the <code>xmldoc</code> 
subdirectory. 
 
<p/> 
In all of the files scanned, lines beginning with 
<code>&lt;flag</code>, <code>&lt;mode</code>, <code>&lt;parm</code>, 
<code>&lt;word</code>, <code>&lt;fvec</code>, <code>&lt;mvec</code>, 
<code>&lt;pvec</code> or <code>&lt;wvec</code> are identified, and the 
information on such a line is used to define a new flag, mode, parameter, 
word, or vector of flags, modes or parameters. To exemplify, consider a line 
<pre> 
&lt;parm name="TimeShower:pTmin" default="0.5" min="0.1" max="2.0"> 
</pre> 
which appears in the <code>TimeShower.xml</code> file, and there 
defines a parameter <code>TimeShower:pTmin</code> with default value 
0.5 GeV and allowed variation in the range 0.1 - 2.0 GeV. The min 
and max values are optional. 
<br/><b>Important:</b> the values in the <code>.xml</code> files should 
not be changed, except by the PYTHIA authors. Any changes should be 
done with the help of the methods described below. 
</li> 
 
<p/> <li> 
Between the creation of the <code>Pythia</code> object and the 
<code>init</code> call for it, you may use several alternative 
methods to modify some of the default values. The same variable 
can be changed several times. If so, it is the last read value 
that counts. The two special 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='Tunes.php?filepath=".$filepath."' target='page'>";?>Tune:ee</a></code> and 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='Tunes.php?filepath=".$filepath."' target='page'>";?>Tune:pp</a></code> modes and the 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='MainProgramSettings.php?filepath=".$filepath."' target='page'>";?>Print:quiet</a></code> flag 
are expanded to change several settings in one go, but these obey 
the same ordering rules. 
 
<p/> 
a) Inside your main program you can directly set values with 
<pre> 
    pythia.readString(string) 
</pre> 
where both the variable name and the value are contained inside 
the character string, separated by blanks and/or a =, e.g. 
<pre> 
    pythia.readString("TimeShower:pTmin = 1.0"); 
</pre> 
The match of the name to the database is case-insensitive. Names 
that do not match an existing variable are ignored. A warning is 
printed, however. Strings beginning with a non-alphanumeric character, 
like # or !, are assumed to be comments and are not processed at all. 
Values below the minimum or above the maximum are set at 
the respective border. In extreme cases, where it is necessary to 
go outside the allowed range, "<code>FORCE=</code>" can replace 
the normal "<code>=</code>" separator to force the requested value, 
at own responsibility. For <code>bool</code> values, the following 
notation may be used interchangeably: 
<code>true = on = yes = ok = 1</code>, while everything else gives 
<code>false</code> (including but not limited to 
<code>false</code>, <code>off</code>, <code>no</code> and 0).<br/> 
 
<p/> 
b) The <code>Pythia</code> <code>readString(string)</code> method 
actually does not do changes itself, but sends on the string either 
to the <code>Settings</code> class or to <code>ParticleData</code>. 
The former holds if the string begins with a letter, the latter 
if it begins with a digit. (The exception is if an input list has been 
begun by an open brace { but no matching close brace } was present; 
then all subsequent non-empty input is directed to <code>Settings</code> 
until the close brace is found.) If desired, it is possible to communicate 
directly with the corresponding <code>Settings</code> method: 
<pre> 
    pythia.settings.readString("TimeShower:pTmin = 1.0"); 
</pre> 
In this case, changes intended for <code>ParticleData</code> 
would not be understood. 
 
<p/> 
c) Underlying the <code>settings.readString(string)</code> method are 
the settings-type-sensitive commands in the <code>Settings</code>, that 
are split by names containing <code>flag</code>, <code>mode</code>, 
<code>parm</code> or <code>word</code>. Thus, the example now reads 
<pre> 
    pythia.settings.parm("TimeShower:pTmin", 1.0); 
</pre> 
Such a form could be convenient e.g. if a parameter is calculated 
at the beginning of the main program, and thus is available as a 
variable rather than as a character string. 
Note that Boolean values must here be given as <code>true</code> or 
<code>false</code> i.e. there is less flexibility than with the 
previous methods. 
 
<p/> 
At the same level, there are several different methods available. 
These are included in the full description below, but normally the user 
should have no need for them. 
 
<p/> 
d) A simpler and more useful way is to collect all your changes 
in a separate file, with one line per change, e.g. 
<pre> 
    TimeShower:pTmin = 1.0 
</pre> 
Each line is read in as a string and processed with the methods already 
introduced. 
 
The file can be read by the 
<pre> 
    pythia.readFile(fileName); 
</pre> 
method (or an <code>istream</code> instead of a <code>fileName</code>). 
The file can freely mix commands to the <code>Settings</code> and 
<code>ParticleData</code> classes, and so is preferable. Lines with 
settings are handled by calls to the 
<code>pythia.settings.readString(string)</code> method. 
 
<p/> 
A file can make use of two extra features that are not available with the 
<code>readString(...)</code> method. One is the possibility to provide 
information for several distinct <?php $filepath = $_GET["filepath"];
echo "<a href='ProgramFlow.php?filepath=".$filepath."' target='page'>";?>subruns</a>. 
The other is the possibility to comment out a section of lines in the file. 
The first line of the commented section should then begin by <code>/*</code> 
and the last begin by <code>*/</code>. This is reminiscent of the convention 
used in C++ and other languages, but is not as powerful, in that it is not 
possible to comment in or out parts of lines. It is only the first two 
non-blank characters of a line that are checked for a match, and a line 
beginning with <code>*/</code> is counted as part of the commented section. 
To avoid mistakes it is best to keep <code>/*</code> and <code>*/</code> 
on lines of their own, optionally followed by comments, but not by commands. 
</li> 
 
<p/> <li> 
In the <code>pythia.init()</code> call, many of the various other program 
elements are initialized, making use of the current values in the database. 
Once initialized, the common <code>Settings</code> database is likely not 
consulted again by these routines. It is therefore not productive to do 
further changes in mid-run: at best nothing changes, at worst you may 
set up inconsistencies. 
 
<p/> 
A routine <code>reInit(fileName)</code> is provided, and can be used to 
zero all the maps and reinitialize them from scratch. Such a call might be 
useful if several subruns are to be made with widely different parameter 
sets - normally the maps are only built from scratch once, namely when the 
<code>Pythia()</code> object is created. A  economical alternative is 
offered by <code>resetAll()</code>, however, which sets all variables back 
to their default values. 
</li> 
 
<p/> <li> 
You may at any time obtain a listing of all variables in the 
database by calling 
<pre> 
    pythia.settings.listAll(); 
</pre> 
The listing is strictly alphabetical, which at least means that names 
from the same file are kept together, but otherwise may not be so 
well-structured: important and unimportant ones will appear mixed. 
A more relevant alternative is 
<pre> 
    pythia.settings.listChanged(); 
</pre> 
where you will only get those variables that differ from their 
defaults. Or you can use 
<pre> 
    pythia.settings.list("string"); 
</pre> 
where only those variables with names that contain the string 
(case-insensitive match) are listed. Thus, with a string 
<code>shower</code>, the shower-related variables would be shown. 
 
<p/> 
The method <code>pythia.settings.output(key)</code> can return the 
value of a variable as a string, convenient for output. In a 
<code>readString</code> or <code>readFile</code> command, the 
construction <code>key = ?</code> will echo back the variable 
and its value, using this method. 
</li> 
 
<p/> <li> 
The above listings are in a tabular form that cannot be read 
back in. Assuming you want to save all changed settings (maybe because 
you read in changes from several files), you can do that by calling 
<pre> 
    pythia.settings.writeFile(fileName); 
</pre> 
This file could then directly be read in by 
<code>readFile(fileName)</code> in a subsequent (identical) run. 
Some variants of this command are listed below. 
</li> 
</ol> 
 
<br/><hr/> 
<a name="section2"></a> 
<h3>The public methods</h3> 
 
The complete list of methods and arguments is as follows. Most of the 
ones of interest to the user have already been mentioned above. 
Others can be used, but the same functionality is better achieved 
by higher-level routines. Some are part of the internal machinery, 
and should not be touched by user. 
 
<p/> 
Note that there is no <code>Settings::readFile(...)</code> method. 
The intention is that you should use <code>Pythia::readFile(...)</code>. 
It parses and decides which individual lines should be sent on to 
<code>Settings::readString(...)</code>. 
 
<a name="anchor1"></a>
<p/><strong> Settings::Settings() &nbsp;</strong> <br/>
the constructor, which takes no arguments. Internal. 
   
 
<a name="anchor2"></a>
<p/><strong> bool Settings::initPtr(Info* infoPtrIn) &nbsp;</strong> <br/>
initialize pointer to error-message database. Internal. 
   
 
<a name="anchor3"></a>
<p/><strong> bool Settings::init(string startFile = &quot;../share/Pythia8/xmldoc/Index.xml&quot;, bool append = false) &nbsp;</strong> <br/>
read in the settings database. 
<br/><code>argument</code><strong> startFile </strong>  : <argument name="startFile" 
default="&quot;../share/Pythia8/xmldoc/Index.xml&quot;"> 
read in the settings from all the files listed in this file, and 
assumed to be located in the same subdirectory. 
   
<br/><code>argument</code><strong> append </strong> (<code>default = <strong>off</strong></code>) :  
By default nothing is done if the method has already been called once. 
If true the further settings read in are added to the current database. 
   
<br/><b>Note:</b> The method returns false if it fails. 
   
 
<a name="anchor4"></a>
<p/><strong> bool init(istream& is, bool append = false) &nbsp;</strong> <br/>
read in the settings from an input stream. This allows initialization 
without reading the xml files directly, which is useful for initialization 
of multiple copies of Pythia8. 
   
 
<a name="anchor5"></a>
<p/><strong> bool Settings::reInit(string startFile = &quot;../share/Pythia8/xmldoc/Index.xml&quot;) &nbsp;</strong> <br/>
overwrite the existing database. 
<br/><code>argument</code><strong> startFile </strong>  : <argument name="startFile" 
default="&quot;../share/Pythia8/xmldoc/Index.xml&quot;"> 
read in the settings from all the files listed in this file, and 
assumed to be located in the same subdirectory. 
   
<br/><b>Note:</b> The method returns false if it fails. 
   
 
<a name="anchor6"></a>
<p/><strong> bool Settings::readString(string line, bool warn = true) &nbsp;</strong> <br/>
read in a string, and change the relevant quantity in the database. 
It is normally used indirectly, via 
<code>Pythia::readString(...)</code> and 
<code>Pythia::readFile(...)</code>. 
<br/><code>argument</code><strong> line </strong>  :  
the string to be interpreted as an instruction. 
   
<br/><code>argument</code><strong> warn </strong> (<code>default = <strong>on</strong></code>) :  
write a warning message or not whenever the instruction does not make 
sense, e.g. if the variable does not exist in the databases. 
   
<br/><b>Note:</b> the method returns false if it fails to 
make sense out of the input string. 
   
 
<a name="anchor7"></a>
<p/><strong> bool Settings::writeFile(string toFile, bool writeAll = false) &nbsp;</strong> <br/>
   
<a name="anchor8"></a>
<strong> bool Settings::writeFile(ostream& os = cout, bool writeAll = false) &nbsp;</strong> <br/>
write current settings to a file or to an <code>ostream</code>. 
<br/><code>argument</code><strong> toFile, os </strong>  :  
file or stream on which settings are written. 
   
<br/><code>argument</code><strong> writeAll </strong> (<code>default = <strong>off</strong></code>) :  
normally only settings that have been changed are written, 
but if true then all settings are output. 
   
<br/><b>Note:</b> the method returns false if it fails. 
   
 
<a name="anchor9"></a>
<p/><strong> bool Settings::writeFileXML(ostream& os = cout) &nbsp;</strong> <br/>
write out the information stored in xmldoc to be used later to 
initialize Settings through an input stream. 
   
 
<a name="anchor10"></a>
<p/><strong> void Settings::listAll() &nbsp;</strong> <br/>
   
<a name="anchor11"></a>
<strong> void Settings::listChanged() &nbsp;</strong> <br/>
   
<a name="anchor12"></a>
<strong> void Settings::list(string match) &nbsp;</strong> <br/>
list all or changed settings, or a group of them. 
<br/><code>argument</code><strong> match </strong>  :  
list all those settings where the name contains 
the <code>match</code> (sub)string (case-insensitive). 
   
   
 
<a name="anchor13"></a>
<p/><strong> string Settings::output(string key, bool fullLine = true) &nbsp;</strong> <br/>
provide the value of a variable as a character string, whatever the type. 
If the variable does not exist then <code>unknown</code> is returned. 
<br/><code>argument</code><strong> key </strong>  :  
the name of the settings variable. 
   
<br/><code>argument</code><strong> fullLine </strong> (<code>default = <strong>on</strong></code>) :  
If true then a whole "line" is returned, " <code>key = value\n</code>", 
while if false only the  <code>value</code> string. 
   
   
 
<a name="anchor14"></a>
<p/><strong> vector&lt;string&gt; Settings::getReadHistory() &nbsp;</strong> <br/>
Method to retrieve the history of <code>readString</code> commands that 
have been processed by the <code>Settings</code> instance, e.g. for 
inspection. Note that <code>readFile</code> command lines are interpreted 
by <code>readString</code> and thus also are listed, as are the 
<code>Settings</code> commands read by <code>Pythia::readString</code> 
and <code>Pythia::readFile</code>. 
   
 
<a name="anchor15"></a>
<p/><strong> vector&lt;string&gt; Settings::getReadHistory(int subrun) &nbsp;</strong> <br/>
Method to retrieve the history of <code>readString</code> commands that 
have been processed by <code>Settings</code>, for a specific subrun 
(see the section on <a href='MainProgramSettings.php?filepath=".$filepath."' target='page'>Main-Program 
Settings</a>). For <code>subrun = -1</code>, returns the 
<code>readString</code> history common to all subruns. For 
<code>subrun >= 0</code>, returns the history of <code>readString</code> 
commands for that specific subrun (omitting the common part). 
   
 
<a name="anchor16"></a>
<p/><strong> void Settings::resetAll() &nbsp;</strong> <br/>
reset all current values to their defaults. 
   
 
<a name="anchor17"></a>
<p/><strong> bool Settings::isFlag(string key) &nbsp;</strong> <br/>
   
<a name="anchor18"></a>
<strong> bool Settings::isMode(string key) &nbsp;</strong> <br/>
   
<a name="anchor19"></a>
<strong> bool Settings::isParm(string key) &nbsp;</strong> <br/>
   
<a name="anchor20"></a>
<strong> bool Settings::isWord(string key) &nbsp;</strong> <br/>
   
<a name="anchor21"></a>
<strong> bool Settings::isFVec(string key) &nbsp;</strong> <br/>
   
<a name="anchor22"></a>
<strong> bool Settings::isMVec(string key) &nbsp;</strong> <br/>
   
<a name="anchor23"></a>
<strong> bool Settings::isPVec(string key) &nbsp;</strong> <br/>
   
<a name="anchor24"></a>
<strong> bool Settings::isWVec(string key) &nbsp;</strong> <br/>
return true if an entry of the given name and kind 
exists, else false. 
   
 
<a name="anchor25"></a>
<p/><strong> void Settings::addFlag(string key, bool default) &nbsp;</strong> <br/>
   
<a name="anchor26"></a>
<strong> void Settings::addMode(string key, int default, bool hasMin, bool hasMax, int min, int max) &nbsp;</strong> <br/>
   
<a name="anchor27"></a>
<strong> void Settings::addParm(string key, double default, bool hasMin, bool hasMax, double min, double max) &nbsp;</strong> <br/>
   
<a name="anchor28"></a>
<strong> void Settings::addWord(string key, string default) &nbsp;</strong> <br/>
   
<a name="anchor29"></a>
<strong> void Settings::addFVec(string key, vector&lt;bool&gt; default) &nbsp;</strong> <br/>
   
<a name="anchor30"></a>
<strong> void Settings::addMVec(string key, vector&lt;int&gt; default, bool hasMin, bool hasMax, int min, int max) &nbsp;</strong> <br/>
   
<a name="anchor31"></a>
<strong> void Settings::addPVec(string key, vector&lt;double&gt; default, bool hasMin, bool hasMax, double min, double max) &nbsp;</strong> <br/>
   
<a name="anchor32"></a>
<strong> void Settings::addWVec(string key, vector&lt;string&gt; default) &nbsp;</strong> <br/>
add an entry of the respective kind to the database. The name and default 
value(s) always has to be supplied, for <code>Mode</code>, <code>Parm</code>, 
<code>MVec</code> and <code>PVec</code> additionally if lower and/or 
upper limits are to be imposed and, if so, what those limit are. 
   
 
<a name="anchor33"></a>
<p/><strong> bool Settings::flag(string key) &nbsp;</strong> <br/>
   
<a name="anchor34"></a>
<strong> int Settings::mode(string key) &nbsp;</strong> <br/>
   
<a name="anchor35"></a>
<strong> double Settings::parm(string key) &nbsp;</strong> <br/>
   
<a name="anchor36"></a>
<strong> string Settings::word(string key) &nbsp;</strong> <br/>
   
<a name="anchor37"></a>
<strong> vector&lt;bool&gt; Settings::fvec(string key) &nbsp;</strong> <br/>
   
<a name="anchor38"></a>
<strong> vector&lt;int&gt; Settings::mvec(string key) &nbsp;</strong> <br/>
   
<a name="anchor39"></a>
<strong> vector&lt;double&gt; Settings::pvec(string key) &nbsp;</strong> <br/>
   
<a name="anchor40"></a>
<strong> vector&lt;string&gt; Settings::wvec(string key) &nbsp;</strong> <br/>
return the current value(s) of the respective setting. If the name 
does not exist in the database, a value <code>false</code>, 
<code>0</code>, <code>0.</code>, <code>&quot; &quot;</code>, or a 
vector of length 1 and value <code>false</code>, <code>0</code>, 
<code>0.</code> or <code>&quot; &quot;</code>, respectively, is returned. 
   
 
<a name="anchor41"></a>
<p/><strong> bool Settings::flagDefault(string key) &nbsp;</strong> <br/>
   
<a name="anchor42"></a>
<strong> int Settings::modeDefault(string key) &nbsp;</strong> <br/>
   
<a name="anchor43"></a>
<strong> double Settings::parmDefault(string key) &nbsp;</strong> <br/>
   
<a name="anchor44"></a>
<strong> string Settings::wordDefault(string key) &nbsp;</strong> <br/>
   
<a name="anchor45"></a>
<strong> vector&lt;bool&gt; Settings::fvecDefault(string key) &nbsp;</strong> <br/>
   
<a name="anchor46"></a>
<strong> vector&lt;int&gt; Settings::mvecDefault(string key) &nbsp;</strong> <br/>
   
<a name="anchor47"></a>
<strong> vector&lt;double&gt; Settings::pvecDefault(string key) &nbsp;</strong> <br/>
   
<a name="anchor48"></a>
<strong> vector&lt;string&gt; Settings::wvecDefault(string key) &nbsp;</strong> <br/>
return the default value(s) of the respective setting. If the name 
does not exist in the database, a value <code>false</code>, 
<code>0</code>, <code>0.</code>, <code>&quot; &quot;</code>, or a 
vector of length 1 and value <code>false</code>, <code>0</code>, 
<code>0.</code> or <code>&quot; &quot;</code>, respectively, is returned. 
   
 
<a name="anchor49"></a>
<p/><strong> map&lt;string, Flag&gt; Settings::getFlagMap(string match) &nbsp;</strong> <br/>
   
<a name="anchor50"></a>
<strong> map&lt;string, Mode&gt; Settings::getModeMap(string match) &nbsp;</strong> <br/>
   
<a name="anchor51"></a>
<strong> map&lt;string, Parm&gt; Settings::getParmMap(string match) &nbsp;</strong> <br/>
   
<a name="anchor52"></a>
<strong> map&lt;string, Word&gt; Settings::getWordMap(string match) &nbsp;</strong> <br/>
   
<a name="anchor53"></a>
<strong> map&lt;string, FVec&gt; Settings::getFVecMap(string match) &nbsp;</strong> <br/>
   
<a name="anchor54"></a>
<strong> map&lt;string, MVec&gt; Settings::getMVecMap(string match) &nbsp;</strong> <br/>
   
<a name="anchor55"></a>
<strong> map&lt;string, PVec&gt; Settings::getPVecMap(string match) &nbsp;</strong> <br/>
   
<a name="anchor56"></a>
<strong> map&lt;string, WVec&gt; Settings::getWVecMap(string match) &nbsp;</strong> <br/>
return a map of all settings of the respective type that contain the 
string "match" in its name. 
   
 
<a name="anchor57"></a>
<p/><strong> void Settings::flag(string key, bool now, bool force = false) &nbsp;</strong> <br/>
   
<a name="anchor58"></a>
<strong> void Settings::mode(string key, int now, bool force = false) &nbsp;</strong> <br/>
   
<a name="anchor59"></a>
<strong> void Settings::parm(string key, double now, bool force = false) &nbsp;</strong> <br/>
   
<a name="anchor60"></a>
<strong> void Settings::word(string key, string now, bool force = false) &nbsp;</strong> <br/>
   
<a name="anchor61"></a>
<strong> void Settings::fvec(string key, vector&lt;bool&gt; now, bool force = false) &nbsp;</strong> <br/>
   
<a name="anchor62"></a>
<strong> void Settings::mvec(string key, vector&lt;int&gt; now, bool force = false) &nbsp;</strong> <br/>
   
<a name="anchor63"></a>
<strong> void Settings::pvec(string key, vector&lt;double&gt; now, bool force = false) &nbsp;</strong> <br/>
   
<a name="anchor64"></a>
<strong> void Settings::wvec(string key, vector&lt;string&gt; now, bool force = false) &nbsp;</strong> <br/>
change the current value(s) of the respective setting to the provided 
new value(s). If lower or upper limits have been set and 
<code>force=false</code>, input values 
outside the allowed range are reinterpreted as being at the nearest 
limit. If <code>force = true</code>, upper and lower limits will be 
ignored, allowing to force values outside the allowed range (to be 
used with caution and at own responsibility!). Any key not found in the 
settings database will be ignored, unless <code>force = true</code>, in 
which case the missing key will be added to the database with the 
given value. 
   
 
<a name="anchor65"></a>
<p/><strong> void Settings::forceMode(string key, int now) &nbsp;</strong> <br/>
   
<a name="anchor66"></a>
<strong> void Settings::forceParm(string key, double now) &nbsp;</strong> <br/>
   
<a name="anchor67"></a>
<strong> void Settings::forceMVec(string key, vector&lt;int&gt; now) &nbsp;</strong> <br/>
   
<a name="anchor68"></a>
<strong> void Settings::forcePVec(string key, vector&lt;double&gt; now) &nbsp;</strong> <br/>
as above, but do not check lower and upper limits, so that the current 
value(s) can be put outside the intended borders. 
<br/><b>Note:</b> these methods have been superseded by the 
<code>force = true</code> option in the standard methods above. 
They are kept for backwards compatibility with version 8.223 and earlier 
but will be removed in a future major release. 
   
 
<a name="anchor69"></a>
<p/><strong> void Settings::resetFlag(string key) &nbsp;</strong> <br/>
   
<a name="anchor70"></a>
<strong> void Settings::resetMode(string key) &nbsp;</strong> <br/>
   
<a name="anchor71"></a>
<strong> void Settings::resetParm(string key) &nbsp;</strong> <br/>
   
<a name="anchor72"></a>
<strong> void Settings::resetWord(string key) &nbsp;</strong> <br/>
   
<a name="anchor73"></a>
<strong> void Settings::resetFVec(string key) &nbsp;</strong> <br/>
   
<a name="anchor74"></a>
<strong> void Settings::resetMVec(string key) &nbsp;</strong> <br/>
   
<a name="anchor75"></a>
<strong> void Settings::resetPVec(string key) &nbsp;</strong> <br/>
   
<a name="anchor76"></a>
<strong> void Settings::resetWVec(string key) &nbsp;</strong> <br/>
reset the current value to the default one. 
   
 
<a name="anchor77"></a>
<p/><strong> bool Settings::getIsInit() &nbsp;</strong> <br/>
return true if the database has been initialized, else false. 
   
 
<a name="anchor78"></a>
<p/><strong> bool Settings::readingFailed() &nbsp;</strong> <br/>
return true if some input could not be parsed, else false. 
   
 
<a name="anchor79"></a>
<p/><strong> bool Settings::unfinishedInput() &nbsp;</strong> <br/>
return true if input of a vector has been begun with am 
open brace { but no matching closing brace } has been found 
(so far), else false. 
   
 
<a name="anchor80"></a>
<p/><strong> bool Settings::onlySoftQCD() &nbsp;</strong> <br/>
return true if the only processes switched on for generation belong 
to the <code>SoftQCD</code> category, else false. Since the method 
contains a hardcoding of what other process types exist, in order to 
detect if any is on, it could be broken by the addition of new internal 
processes. It is therefore mainly useful for warning purposes, not for 
hard decisions. 
   
 
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
