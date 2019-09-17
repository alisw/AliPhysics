<html>
<head>
<title>Histograms</title>
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

<form method='post' action='Histograms.php'>
 
<h2>Histograms</h2> 
<ol id="toc">
  <li><a href="#section0">Basic principles</a></li>
  <li><a href="#section1">Basic output format</a></li>
  <li><a href="#section2">Matplotlib output format</a></li>
  <li><a href="#section3">The methods</a></li>
</ol>

 
The <code>Hist</code> class gives a simple implementation of 
one-dimensional histograms, useful for quick-and-dirty testing, 
without the need to link to more sophisticated packages. 
For this reason it is used in many of the 
<?php $filepath = $_GET["filepath"];
echo "<a href='SampleMainPrograms.php?filepath=".$filepath."' target='page'>";?>sample main programs</a> 
found in the <code>examples</code> subdirectory. 
 
<a name="section0"></a> 
<h3>Basic principles</h3> 
 
We here provide a simple overview of what is involved. 
As a first step you need to declare a histogram, with name, 
title, number of bins and <i>x</i> range (from, to). 
<pre> 
   Hist ZpT( "Z0 pT spectrum", 100, 0., 100.); 
</pre> 
Alternatively you can first declare it and later define it: 
<pre> 
   Hist ZpT; 
   ZpT.book( "Z0 pT spectrum", 100, 0., 100.); 
</pre> 
 
Once declared, its contents can be added by repeated calls to 
<code>fill</code>, 
<pre> 
   ZpT.fill( 22.7, 1.); 
</pre> 
where the first argument is the <i>x</i> value and the second the 
weight. Since the weight defaults to 1 the last argument could have 
been omitted in this case. 
 
<p/> 
A set of overloaded operators have been defined, so that histograms 
can be added, subtracted, divided or multiplied by each other. Then the 
contents are modified accordingly bin by bin. Thus the relative 
deviation between two histograms <code>data</code> and 
<code>theory</code> can be found as 
<pre> 
  diff = (data - theory) / (data + theory); 
</pre> 
assuming that <code>diff</code>, <code>data</code> and <code>theory</code> 
have been booked with the same number of bins and <i>x</i> range. That 
responsibility rests on the user; some checks are made for compatibility, 
but not enough to catch all possible mistakes. 
 
<p/> 
Also overloaded operations with double real numbers are available. 
Again these four operations are defined bin by bin, i.e. the 
corresponding amount is added to, subtracted from, multiplied by or 
divided by each bin. The double number can come before or after the 
histograms, with obvious results. Thus the inverse of a histogram 
<code>result</code> is given by <code>1. / result</code>. 
The two kind of operations can be combined, e.g. 
<pre> 
  allpT = ZpT + 2. * WpT 
</pre> 
Finally, also the <code>+=, -+, *=, /=</code> are overloaded, with 
the right-hand side being either a histogram or a real number. 
 
<a name="section1"></a> 
<h3>Basic output format</h3> 
 
A histogram can be printed by making use of the overloaded &lt;&lt; 
operator, e.g.: 
<pre> 
   cout &lt;&lt; ZpT; 
</pre> 
The printout format is inspired by the old HBOOK one. To understand 
how to read this format, consider the simplified example 
<pre> 
        3.50*10^ 2  9 
        3.00*10^ 2  X   7 
        2.50*10^ 2  X  1X 
        2.00*10^ 2  X6 XX 
        1.50*10^ 2  XX5XX 
        1.00*10^ 2  XXXXX 
        0.50*10^ 2  XXXXX 
 
          Contents 
            *10^ 2  31122 
            *10^ 1  47208 
            *10^ 0  79373 
 
          Low edge  -- 
            *10^ 1  10001 
            *10^ 0  05050 
</pre> 
The key feature is that the <code>Contents</code> and 
<code>Low edge</code> have to be read vertically. For instance, 
the first bin has the contents 
<code>3 * 10^2 + 4 * 10^1 + 7 * 10^0 = 347</code>. Correspondingly, 
the other bins have contents 179, 123, 207 and 283. The first bin 
stretches from <code>-(1 * 10^1 + 0 * 10^0) = -10</code> to the 
beginning of the second bin, at <code>-(0 * 10^1 + 5 * 10^0) = -5</code>. 
 
<p/> 
The visual representation above the contents give a simple impression 
of the shape. An <code>X</code> means that the contents are filled up 
to this level, a digit in the topmost row the fraction to which the 
last level is filled. So the 9 of the first column indicates this bin 
is filled 9/10 of the way from <code>3.00*10^2 = 300</code> to 
<code>3.50*10^2 = 350</code>, i.e. somewhere close to 345, 
or more precisely in the range 342.5 to 347.5. 
 
<p/> 
The printout also provides some other information, such as the 
number of entries, i.e. how many times the histogram has been filled, 
the total weight inside the histogram, the total weight in underflow 
and overflow, and the mean value and root-mean-square width (disregarding 
underflow and overflow). The mean and width assumes that all the 
contents is in the middle of the respective bin. This is especially 
relevant when you plot a integer quantity, such as a multiplicity. 
Then it makes sense to book with limits that are half-integers, e.g. 
<pre> 
   Hist multMPI( "number of multiparton interactions", 20, -0.5, 19.5); 
</pre> 
so that the bins are centered at 0, 1, 2, ..., respectively.  This also 
avoids ambiguities which bin gets to be filled if entries are 
exactly at the border between two bins. Also note that the 
<code>fill( xValue)</code> method automatically performs a cast 
to double precision where necessary, i.e. <code>xValue</code> 
can be an integer. 
 
<a name="section2"></a> 
<h3>Matplotlib output format</h3> 
 
Assuming you have Python installed on your platform, it is possible to 
generate simple <a href="https://matplotlib.org/">Matplotlib/Pyplot</a> 
Python code from the histograms generated above, which then can be run 
to produce PDF plots. This should be done near the end of a run, after 
the histograms have been filled and properly normalized, as an alternative 
or complement to the basic output format above. 
 
<p/> 
In  a first step you must then decide on the name of the Python program, 
e.g.: 
<pre> 
  HistPlot hpl( "bosonpT"); 
</pre> 
where file ending <code>.py</code> is added automatically. 
 
<p/> 
For each new frame the name should be given, which will later give 
rise to a PDF file, with ending <code>.pdf</code> added automatically. 
If you leave the name field empty the same file will be used as for the 
latest named one, i.e. producing several frames in one file. One can 
optionally also give title and <i>x</i> and <i>y</i> axis labels: 
<pre> 
  hpl.frame( "pTdist", "Boson pT distributions", "pT (GeV)", "sigma"); 
</pre> 
Next, existing <code>Hist</code> histograms can be added to the frame, 
one by one: 
<pre> 
  hpl.add( ZpT, "-"); 
  hpl.add( WpT, "--,indigo"); 
  hpl.add( ZpT + 2. * WpT , "", "pT spectrum of Z, W+ and W-"); 
</pre> 
where the second argument tells how each histogram will be plotted. 
Default is histogram style, "h", but the values can also be connected 
with full lines "-", dashed ones "--", or dash-dotted ones "-.", or 
plotted as points "." or crosses "x", to mention some of the many 
options offered by Pyplot. Here you can also specify the colour, 
separated by a comma from the line style, to override the normal 
colour cycle. The most common colours can be given just as a single 
letter, such as "r", "g", "b", but a more extensive 
<a href="https://matplotlib.org/examples/color/named_colors.html"> 
colour palette</a> allow finetuning to nuances such as "orange", 
"gold", "darkgreen", "royalblue", "crimson", and so on. 
A third argument can set the legend of each histogram; by default 
it is taken as the title of histogram. 
<br/>Finally the plot code itself will be set up by 
<pre> 
  hpl.plot(); 
</pre> 
where optionally it is possible to demand a logarithmic <i>y</i> scale. 
 
<p/> 
The <code>frame - add - plot</code> steps can be repeated as needed, each 
giving rise to a separate PDF file with a plot. In case a plot is to be 
generated from a single histogram the three steps can be joined into one 
<pre> 
  hpl.plotFrame( "onlyZ", ZpT ); 
</pre> 
where only the name of the PDF file and the histogram are compulsory, while 
further arguments as discussed above are optional. 
 
<p/> 
At the end, a file <code>bosonpT.py</code> has now been generated with the 
proper plotting commands. Additionally a data file has been generated for 
each histogram to be plotted, <code>pTdist-0.dat, pTdist-1.dat</code>, etc. 
Now doing 
<pre> 
  python bosonpT.py 
</pre> 
in a terminal window will produce the plots, such as <code>pTdist.pdf</code> 
and <code>onlyZ.pdf</code>. You may of course edit the python file further 
to improve on the finer details. 
 
<p/> 
Examples are provided in <code>main03.cc</code> and <code>main07.cc</code>, 
where the latter is the simpler one. The <code>main51.cc</code> example 
illustrates that the <i>x</i> scale can be chosen logarithmically, 
by using an optional last argument when histograms are booked. 
 
<a name="section3"></a> 
<h3>The methods</h3> 
 
We here collect a more complete and formal overview of the methods. 
 
<a name="anchor1"></a>
<p/><strong> Hist::Hist() &nbsp;</strong> <br/>
declare a histogram, but does not define it. 
   
 
<a name="anchor2"></a>
<p/><strong> Hist::Hist(string title, int numberOfBins, double xMin, double xMax, bool logX = false) &nbsp;</strong> <br/>
declare and define a histogram, where 
<br/><code>argument</code><strong> title </strong>  :  
is a string with the title of the histogram at output, 
   
<br/><code>argument</code><strong> numberOfBins </strong>  :  
is the number of bin the <i>x</i> range will be subdivided into, 
limited to be at most 1000, 
   
<br/><code>argument</code><strong> xMin </strong>  :  
is the lower edge of the histogram, 
   
<br/><code>argument</code><strong> xMax </strong>  :  
is the upper edge of the histogram, 
   
<br/><code>argument</code><strong> logX </strong>  :  
gives a logarithmic <i>x</i> scale between <code>xMin</code> and 
<code>xMax</code> if set <code>true</code>. Note that then 
<code>xMin</code> must be above 1e-20 so as to stay with positive 
numbers. 
   
   
 
<a name="anchor3"></a>
<p/><strong> Hist::Hist(const Hist&amp; h) &nbsp;</strong> <br/>
creates an identical copy of the histogram in the argument, 
including bin contents. 
   
 
<a name="anchor4"></a>
<p/><strong> Hist::Hist(string title, const Hist&amp; h) &nbsp;</strong> <br/>
creates an identical copy of the histogram in the argument, 
including bin contents, except that a new title is provided 
as first argument. 
   
 
<a name="anchor5"></a>
<p/><strong> Hist&amp; Hist::operator=(const Hist&amp; h) &nbsp;</strong> <br/>
copies all properties of the histogram in the argument, 
except that the original histogram title is retained. 
   
 
<a name="anchor6"></a>
<p/><strong> void Hist::book(string title, int numberOfBins, double xMin, double xMax, bool logXIn = false) &nbsp;</strong> <br/>
define a histogram that previously was only declared; 
see above for the meaning of the arguments. 
   
 
<a name="anchor7"></a>
<p/><strong> void Hist::title(string title) &nbsp;</strong> <br/>
change the title of a histogram, but keep other properties unchanged. 
   
 
<a name="anchor8"></a>
<p/><strong> void Hist::null() &nbsp;</strong> <br/>
reset bin contents, but keep other histogram properties unchanged. 
   
 
<a name="anchor9"></a>
<p/><strong> void Hist::fill(double xValue, double weight) &nbsp;</strong> <br/>
fill the histogram, where 
<br/><code>argument</code><strong> xValue </strong>  :  
is the <i>x</i> position where the filling should occur, and 
   
<br/><code>argument</code><strong> weight </strong> (<code>default = <strong>1.</strong></code>) :  
is the amount of weight to be added at this <i>x</i> value. 
   
   
 
<a name="anchor10"></a>
<p/><strong> friend ostream& operator&lt;&lt;(ostream&amp; os, const Hist&amp; h) &nbsp;</strong> <br/>
appends a simple histogram printout (see above for format) to the 
<code>ostream</code>, while leaving the histogram object itself 
unchanged. At most 100 columns are allowed to be displayed. 
If the number of bins is larger than 100 then the contents of 
adjacent bins are added to give the value in each column. (Two by two 
up to 200 bins, three by three up to 300, and so on, with the very 
last column possibly summing fewer rows than the others.) 
<br/>If the histogram has been booked with logarithmic <i>x</i> 
scale then the log10 of each lower bin border will be displayed 
rather than the border itself, to avoid undue complication. 
The Pyplot interface will display correct values, however. 
   
 
<a name="anchor11"></a>
<p/><strong> void Hist::table(ostream&amp; os = cout, bool printOverUnder = false, bool xMidBin = true) &nbsp;</strong> <br/>
   
<a name="anchor12"></a>
<strong> void Hist::table(string fileName, bool printOverUnder = false, bool xMidBin = true) &nbsp;</strong> <br/>
print a two-column table, where the first column gives the center of 
each bin and the second one the corresponding bin contents. The table 
may be useful for plotting e.g. with Gnuplot. 
<br/>The desired output stream or file name can be provided as argument. 
The former is more flexible (e.g., it allows easy append to an existing 
file), whereas the latter is simpler for the case that each histogram 
should be a file of its own. 
<br/>An optional <code>printOverUnder = true</code> argument allows also 
underflow and overflow contents to be printed. (The arbitrary <i>x</i> 
coordinates for these are placed as if corresponding to same-size bins 
just below or above the regular histogram bins.) 
<br/>An optional <code>xMidBin = false</code> argument will have the 
<i>x</i> value at the beginning of each bin printed, rather than the 
default midpoint value. 
   
 
<a name="anchor13"></a>
<p/><strong> void Hist::rivetTable(ostream&amp; os = cout, bool printError = false) &nbsp;</strong> <br/>
   
<a name="anchor14"></a>
<strong> void Hist::rivetTable(string fileName, bool printError = false) &nbsp;</strong> <br/>
print a five-column table, where the first two columns give the lower 
and upper borders of each bin, the third one the bin contents, and the 
fourth and fifth the error (up and down) associated with the contents. 
This format matches the one that Rivet uses for its histograms. 
The choice between the two methods is the same as above for the 
<code>table</code> methods. 
<br/>The error bins are put to zero by default, since the PYTHIA 
histogramming is not sophisticated enough to compensate for rescalings 
or other operations, or for weighted events. With the optional 
<code>printError = true</code> the error will be taken as the 
square root of the bin content, as is relevant if this content has 
the same unit weight for each entry to it. 
   
 
<a name="anchor15"></a>
<p/><strong> void Hist::pyplotTable(ostream&amp; os = cout, bool isHist = true) &nbsp;</strong> <br/>
   
<a name="anchor16"></a>
<strong> void Hist:pyplotTable(string fileName, bool isHist = true) &nbsp;</strong> <br/>
prints either a two- or a three-column table, depending on whether 
the contents are to be plotted as a curve or as a histogram, the latter 
being default. In either case the first column gives the center of each 
bin and the second one the corresponding bin contents. For a histogram 
the third column gives the lower border of each bin. Then also a final 
line is provided, that gives no further histogram contents but provides 
the upper border of the last bin. This format is used for plotting with 
Pyplot, see below. The choice between the two methods is the same as 
above for the <code>table</code> methods. 
   
 
<a name="anchor17"></a>
<p/><strong> friend void table(const Hist&amp; h1, const Hist&amp; h2, ostream& os = cout, bool printOverUnder = false, bool xMidBin = true) &nbsp;</strong> <br/>
   
<a name="anchor18"></a>
<strong> friend void table(const Hist&amp; h1, const Hist&amp; h2, string fileName, bool printOverUnder = false, bool xMidBin = true) &nbsp;</strong> <br/>
print a three-column table, where the first column gives the center of 
each bin and the second and third ones the corresponding bin contents 
of the two histograms. Only works if the two histograms have the same 
<i>x</i> axis (within a tiny tolerance), else nothing will be done. 
The optional last two arguments allows also underflow and overflow 
contents to be printed, and the <i>x</i> to refer to the beginning 
of the bin rather than the center; see above. 
   
 
<a name="anchor19"></a>
<p/><strong> string Hist::getTitle() &nbsp;</strong> <br/>
return the title of the histogram. 
   
 
<a name="anchor20"></a>
<p/><strong> int Hist::getBinNumber() &nbsp;</strong> <br/>
return the number of bins in the histogram. 
   
 
<a name="anchor21"></a>
<p/><strong> bool Hist::getLinX() &nbsp;</strong> <br/>
return true if the histogram has a linear <i>x</i> scale and 
false if a logarithmic one. 
   
 
<a name="anchor22"></a>
<p/><strong> double Hist::getBinContent(int iBin) &nbsp;</strong> <br/>
return the value in bin <code>iBin</code>, ranging from 1 through 
<code>numberOfBins</code>, with <code>0</code> for underflow and 
<code>numberOfBins + 1</code> for overflow. 
   
 
<a name="anchor23"></a>
<p/><strong> int Hist::getEntries() &nbsp;</strong> <br/>
return the number of entries, i.e. the number of time that 
<code>fill(...)</code> has been called. 
   
 
<a name="anchor24"></a>
<p/><strong> bool Hist::sameSize(const Hist&amp; h) &nbsp;</strong> <br/>
checks that the number of bins and upper and lower limits are the 
same as in the histogram in the argument. 
   
 
<a name="anchor25"></a>
<p/><strong> void Hist::takeLog(bool tenLog = true) &nbsp;</strong> <br/>
by default take 10-logarithm of current contents bin by bin. With 
optional argument <code>false</code> instead take <i>e</i>-logarithm 
of contents bin by bin. If to be used, then right before the 
histogram is output. 
   
 
<a name="anchor26"></a>
<p/><strong> void Hist::takeSqrt() &nbsp;</strong> <br/>
take square root of current contents bin by bin, with negative contents 
set to zero. 
   
 
<a name="anchor27"></a>
<p/><strong> double Hist::smallestAbsValue() &nbsp;</strong> <br/>
returns the smallest absolute value (above 1e-20, to sidestep bins with zero 
content) in any of the bins. 
   
 
<a name="anchor28"></a>
<p/><strong> Hist&amp; Hist::operator+=(const Hist&amp; h) &nbsp;</strong> <br/>
   
<a name="anchor29"></a>
<strong> Hist&amp; Hist::operator-=(const Hist&amp; h) &nbsp;</strong> <br/>
adds or subtracts the current histogram by the contents of the 
histogram in the argument if <code>sameSize(...)</code> is true, 
else does nothing. 
   
 
<a name="anchor30"></a>
<p/><strong> Hist&amp; Hist::operator*=(const Hist&amp; h) &nbsp;</strong> <br/>
   
<a name="anchor31"></a>
<strong> Hist&amp; Hist::operator/=(const Hist&amp; h) &nbsp;</strong> <br/>
multiplies or divides the current histogram by the contents of the 
histogram in the argument if <code>sameSize(...)</code> is true, 
else does nothing. 
   
 
<a name="anchor32"></a>
<p/><strong> Hist&amp; Hist::operator+=(double f) &nbsp;</strong> <br/>
   
<a name="anchor33"></a>
<strong> Hist&amp; Hist::operator-=(double f) &nbsp;</strong> <br/>
adds or subtracts each bin content by the common offset <i>f</i>. 
   
 
<a name="anchor34"></a>
<p/><strong> Hist&amp; Hist::operator*=(double f) &nbsp;</strong> <br/>
   
<a name="anchor35"></a>
<strong> Hist&amp; Hist::operator*=(double f) &nbsp;</strong> <br/>
multiplies or divides each bin content by the common factor <i>f</i>. 
   
 
<a name="anchor36"></a>
<p/><strong> friend Hist operator+(double f, const Hist&amp; h1) &nbsp;</strong> <br/>
   
<a name="anchor37"></a>
<strong> friend Hist operator+(const Hist&amp; h1, double f) &nbsp;</strong> <br/>
   
<a name="anchor38"></a>
<strong> friend Hist operator+(const Hist&amp; h1, const Hist h2) &nbsp;</strong> <br/>
add a constant to a histogram or two histograms to each other, bin by bin. 
   
 
<a name="anchor39"></a>
<p/><strong> friend Hist operator-(double f, const Hist&amp; h1) &nbsp;</strong> <br/>
   
<a name="anchor40"></a>
<strong> friend Hist operator-(const Hist&amp; h1, double f) &nbsp;</strong> <br/>
   
<a name="anchor41"></a>
<strong> friend Hist operator-(const Hist&amp; h1, const Hist h2) &nbsp;</strong> <br/>
subtract a histogram from a constant, a constant from a histogram, 
or two histograms from each other, bin by bin. 
   
 
<a name="anchor42"></a>
<p/><strong> friend Hist operator*(double f, const Hist&amp; h1) &nbsp;</strong> <br/>
   
<a name="anchor43"></a>
<strong> friend Hist operator*(const Hist&amp; h1, double f) &nbsp;</strong> <br/>
   
<a name="anchor44"></a>
<strong> friend Hist operator*(const Hist&amp; h1, const Hist h2) &nbsp;</strong> <br/>
multiply a constant by a histogram or two histograms by each other, 
bin by bin. 
   
 
<a name="anchor45"></a>
<p/><strong> friend Hist operator/(double f, const Hist&amp; h1) &nbsp;</strong> <br/>
   
<a name="anchor46"></a>
<strong> friend Hist operator/(const Hist&amp; h1, double f) &nbsp;</strong> <br/>
   
<a name="anchor47"></a>
<strong> friend Hist operator/(const Hist&amp; h1, const Hist h2) &nbsp;</strong> <br/>
divide a constant by a histogram, a histogram by a constant, 
or two histograms by each other, bin by bin. 
   
 
<a name="anchor48"></a>
<p/><strong> HistPlot::HistPlot(string pythonName) &nbsp;</strong> <br/>
create a file to which successively Python commands can be written. 
<br/><code>argument</code><strong> pythonName </strong>  :  
the name of the Python code file will become <code>pythonName.py</code>. 
When later executed with <code>python pythonName.py</code> the frames 
defined below will be drawn and output. 
   
   
 
<a name="anchor49"></a>
<p/><strong> void HistPlot::frame( string frameName, string title = &quot;&quot;, string xLabel = &quot;&quot;, string yLabel = &quot;&quot;) &nbsp;</strong> <br/>
command to open a frame where successive histograms can be added and 
plotted, where 
<br/><code>argument</code><strong> frameName </strong>  :  
is the frame that will become a PFD file <code>frameName.pdf</code> when 
the Python program is run; if left blank the frame will be appended to the 
most recent previous  PDF file, so clearly it cannot be left blank for the 
first call, 
   
<br/><code>argument</code><strong> title </strong>  :  
is an optional title of the histogram, plotted on top of the frame, 
   
<br/><code>argument</code><strong> xLabel </strong>  :  
is an optional label plotted below the <i>x</i> axis, and 
   
<br/><code>argument</code><strong> yLabel </strong>  :  
is an optional label plotted to the right of the <i>y</i> axis. 
   
<br/><b>Note:</b> the title and labels may contain special LaTeX symbols, 
e.g. simple math formulae enclosed by dollar signs $...$, that are 
interpreted by the Pyplot package. A backslash \ needs to be doubled to \\, 
however, since else it will be intepreted as an escape sequence. 
   
 
<a name="anchor50"></a>
<p/><strong> void HistPlot::add(const Hist& hist, string style = &quot;h&quot;, string legend = &quot;void&quot;) &nbsp;</strong> <br/>
add a histogram to the current plot frame, where 
<br/><code>argument</code><strong> hist </strong>  :  
is the name of the histogram to be plotted, 
   
<br/><code>argument</code><strong> style </strong>  :  
tells how the histogram data will be represented in the figure, where 
default "h" is histogram style, "-" full lines, "--" dashed ones, 
"-." dash-dotted ones, "." points, "o" circles, "+" and "x" crosses, 
"*" stars, and so on as described in the Pyplot manual; additionally a 
colour may also be specified as described earlier, separated by a comma 
from the symbol; either of both can be used, e.g. "h,r", ",magenta", 
"+,", "--" or just ""; and 
   
<br/><code>argument</code><strong> legend </strong>  :  
is the text written in association with the plotting symbol, where the 
default <code>void</code> is replaced by the title provided when the 
histogram was booked, and LaTex input may be used as described in the note 
above. 
   
   
 
<a name="anchor51"></a>
<p/><strong> void HistPlot::plot( bool logY = false) &nbsp;</strong> <br/>
write the Pyplot commands for the current frame, as defined by the input 
from the two methods above, and additionally write the data files needed 
when the Python code is to be run. 
<br/><code>argument</code><strong> logY </strong>  :  
By default a linear <i>y</i> scale is used, but when true the <i>y</i> 
scale is instead logarithmic, except that a linear scale is used the very 
last step to zero bin content. 
   
<br/><b>Note:</b> whether the <i>x</i> scale is linear or logarithmic 
depends on what choice was made when the histogram was booked, and cannot 
be changed later. 
   
 
<a name="anchor52"></a>
<p/><strong> void HistPlot::plotFrame( string frameName, const Hist& hist, string title = &quot;&quot;, string xLabel = &quot;&quot;, string yLabel = &quot;&quot;, string style = &quot;h&quot;, string legend = &quot;void&quot;, bool logY = false) &nbsp;</strong> <br/>
is an omnibus version combining the three above methods when only a single 
histogram is to be plotted in a frame. The arguments and their meaning is 
the same as already explained for them, and the order in which they appear 
is almost the same, the exception being that the two compulsory arguments 
<code>frameName</code> and <code>hist</code> have been put at the beginning. 
   
 
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
