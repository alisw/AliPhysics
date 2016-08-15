# A class to hold results with statistical and systematic errors {#gse_index}

## The class

The class _GraphSysErr_ lives in `GraphSysErr.C`.  The class is
heavily commented and marked up for documentation using Doxygen.

## Features

- Support multiple types of systematic uncertainties
  - Uncertainties common (correlated) between all points
  - Uncertainties for each individual (uncorrelated) point 
- Support statistical errors.
- Uncertainties can be relative or absolute
- Import data sets from HepData input format
- Export data to HepData input format
  - Sets of objects can be exported as tables
  - Number of significant digits can be specified 
- Visual representation in many ways 
- Rigorous error propagation 
- Ratios, scaling, etc. of graphs 
- Combination (averages) of data sets.
- ... 

## The tests

- TestGSE.C runs an extensive test of the various plotting options.
- TestImport.C Import a single dataset from Durham database
- TestSample.C read the canonical sample input file 

## Installation

Get the current tar-ball from

* http://cern.ch/cholm/root/gse/GSE-current.tar.gz

unpack and copy the file `GraphSysErr.C` anywhere you like.

Alternatively, one can check out the source from the subversion
repository

    svn co https://svn.cern.ch/reps/alicefmd/graphsyserr/trunk graphsyserr

Test scripts and data is available in the tar-ball together with this
information.  To run a full test, do

    make test

This will run the `TestGSE.C` example, and then import the generated
`test.input` using `TestImport.C` and then re-export to `check.input`
and finally compare the two generated files. 

## Usage

First of all, we should load the class.  It is best if we AcLic
compile the code

~~~{.cxx}
    if (!gROOT->GetClass("GraphSysErr",false,true))
      gROOT->LoadMacro("GraphSysErr.C+g");
~~~

### Declare an object

First, one declares an object of the class GraphSysErr with a name and
a title.  Optionally, one can specify how many data points will be
present.

~~~{.cxx}
    GraphSysErr* gse = new GraphSysErr("data","My Data");
~~~
	
### Declare the systematic errors

We then declare which systematic errors we have. Systematic errors are
either _common_ or _point-to-point_. Furthermore, each type of
systematic error can be specified as either relative (to the point
value) or absolute. 

~~~{.cxx}
    // Two sources of common errors one relative, one absolue
    UInt_t cm1 = gse->DefineCommon("Common 0.05", false, .05);
    UInt_t cm2 = gse->DefineCommon("Common 10%", true, .1);
    // Two sources of point-to-point errors, one relative, one absolute
    UInt_t pp1 = gse->DeclarePoint2Point("Point-to-Point 0.1-0.2", true);
    UInt_t pp2 = gse->DeclarePoint2Point("Point-to-Point 5-10%", false)
~~~
	
The unsigned integer returned by these member functions are abstracted
unique (to each object) identifiers.  This identifier is used wherever
one need to manipulate the systematic errors.

### Specify the data set and systematic errors

We can now specify the data points, statistical, and the
point-to-point systematic errors.  For this example we sample a
Gaussian distribution and we take the bin center and frequency as the
point, and the square root of the number of counts as the statistical
error

~~~{.cxx}
    // Fill a histogram with a Guassian random deviate 
    TH1* h = new TH1F("h", "h", 30, -3, 3);
    h->Sumw2();
    h->SetDirectory(0);
    h->FillRandom("gaus",1000);
    h->Scale(1./1000, "width");

    for (Int_t i = 0; i < h->GetNbinsX(); i++) { 
      Int_t    bin = i+1;
      Double_t x   = h->GetXaxis()->GetBinCenter(bin);
      Double_t y   = h->GetBinContent(bin);
      Double_t sta = h->GetBinError(bin);
      
      // Set data 
      gse->SetPoint(i, x, y);
      gse->SetStatError(i, sta);

      // Set point-to-point errors 
      gse->SetSysError(pp1, i, 0., gRandom->Uniform(0.1, 0.2));
      gse->SetSysError(pp2, i, 0., gRandom->Uniform(0.05, 0.1));
    }
~~~

We can now plot, store, or export this data set.  However, we should
customize the plot a little before that.


### Styling the plot

First, the overall options for the data and statistical errors it
self.

~~~{.cxx}
    // Draw data with-out ticks 
    gse->SetDataOption(GraphSysErr::kNoTick);
	// Set axis titles
	gse->SetXTitle("X");
	gse->SetYTitle("Y");
~~~
	
Then, in case we plot the combined errors, set the appearance 

~~~{.cxx}
    // Set options on summed errors (in case of option COMBINED)
    gse->SetSumLineColor(kRed+2);
    gse->SetSumLineWidth(2);
    gse->SetSumTitle("All errors");
    gse->SetSumOption(GraphSysErr::kHat);
~~~
	
Set the style of the common systematic errors:

~~~{.cxx}
    // Set attributes of common errors 
    gse->SetSysFillColor(cm1, kRed+2);
    gse->SetSysFillStyle(cm1, 3001);
    gse->SetSysLineColor(cm1, kRed+2);
    gse->SetSysFillColor(cm2, kCyan+2);
    gse->SetSysFillStyle(cm2, 3001);
    gse->SetSysOption(cm1, GraphSysErr::kBox);
    gse->SetSysOption(cm2, GraphSysErr::kRect);
~~~
	
Set the style on the point-to-point systematic errors

~~~{.cxx}
    // Set attributes of other errors 
    gse->SetSysLineColor(pp1, kBlue+2);
    gse->SetSysLineWidth(pp1, 2);
    gse->SetSysLineColor(pp2, kGreen+2);
    gse->SetSysLineWidth(pp2, 3);
    gse->SetSysOption(pp1, GraphSysErr::kHat);
    gse->SetSysOption(pp2, GraphSysErr::kBar);    
~~~
	
Note that by default, and if one does not specify errors along X, the
horizontal size of the `GraphSysErr::kBox` errors are the same as the
vertical size.  To limit the size, set the `TStyle::SetErrorX` option:

~~~{.cxx}
    // Adjust the horizontal size of boxes and rectangles if no X
	// errors are set. 
    gStyle->SetErrorX(.2);
    // Adjust size of hat, cap, ...
    gStyle->SetEndErrorSize(6);
~~~

### Draw the plot

The data can be plotted in a variety of ways.  We can combine all
systematic point-to-point systematic errors, either as the square root
of the sum of the squares, or as direct sums. 

~~~{.cxx}
    gse->Draw("COMBINE QUADRATURE AXIS MAX");
	gse->Draw("COMBINE DIRECT AXIS MAX");
~~~
	
The statistical errors and the common systematic errors will be
displaced separately. The option MAX indicates that the common
systematic errors will be displaced near the maximum of the data.  The
option AXIS means we will draw the XY axis frame.  With out it, the
data will be drawn in the current axis frame, 

To combine all systematic errors - including common systematic
errors - we add COMMON to the options

~~~{.cxx}
    gse->Draw("COMBINE COMMON QUADRATURE AXIS MAX");
	gse->Draw("COMBINE COMMON DIRECT AXIS MAX");
~~~
	
If we also want to combine the statistical errors, we add the option
STAT 

~~~{.cxx}
    gse->Draw("COMBINE COMMON STAT QUADRATURE AXIS MAX");
	gse->Draw("COMBINE COMMON STAT DIRECT AXIS MAX");
~~~
	
We can also stack up the systematic errors.  To stack up
point-to-point systematic errors but display the common systematic and
statistical errors separately we do 


~~~{.cxx}
    gse->Draw("STACK QUADRATURE AXIS MIN");
	gse->Draw("STACK DIRECT AXIS MIN");
~~~
	
Here, MIN means display the common systematic errors near the minimum
of the data.  Again, with options COMMON and STAT we can also stack up
the other kinds of errors 

~~~{.cxx}
    gse->Draw("STACK COMMON QUADRATURE AXIS MIN");
	gse->Draw("STACK COMMON STAT DIRECT AXIS MIN");
~~~
	
### Storing the result

We can store the object as any other ROOT object by simply writing it
to a ROOT file

~~~{.cxx}
    TFile* file = TFile::Open("output.root", "RECREATE");
	gse->Write();
    file->Close();
~~~
### Export to file for upload to Durham database

A useful feature of this class is that we can export the data to a
_Durham database_ input formatted file (see also BNF.md).  This file
can, after a little bit of editing, be uploaded directly to the Durham
database.  To export a single data set to standard output we do

~~~{.cxx}
    gse->Export(false);
~~~

If the first argument is `false` we only output the data set.  If it
is `true` we also output header information.  To write to a specific
file, first open a stream to that file, and pass that object as the
second argument. 

~~~{.cxx}
    std::ofstream output("data.input");
	gse->Export(true,output);
	output.close();
~~~
	
Suppose one has a collection `allData` with some data sets to
exported, we can simply iterate over that list and write to an output
file.

~~~{.cxx}
	TList* allData = new TList;
    GraphSysErr* first = new GraphSysErr;
    first->SetKey("author", "Last name of first author");
    first->SetKey("reference", "Journal/archive reference");
    first->SetKey("doi", "DOI Id");
    first->SetKey("laboratory", "e.g., CERN");
    first->SetKey("accelerator", e.g., Tevatron");
    first->SetKey("experiment", "e.g., STAR");
    first->SetKey("inspireId", "Reference to inSpire");
    first->SetKey("cdsId", "CERN Document Server Id");
    first->SetKey("title", "Article title");
    first->SetKey("abstract", "Short summery of data");
    allData->Add(first);

    // Add more graphs to allData

    std::ofstream output("article.input");
	TIter next(allData);
	GraphSysErr* g = 0;
	Bool_t first = true;
	while ((g = static_cast<GraphSysError*>(next()))) {
	  g->Export(output, (first ? "HC" : ""));
	  first = false;
    }
~~~


Suppose we have a list `allCuts` of related plots we want to export to
a single table, we can do

~~~{.cxx}
    std::ofstream table("table.input");
    GraphSysErr::Export(allCuts, table);
	table.close();
~~~


And similarly suppose we have a list `allPlots` of all plots in a
paper, and where some of them are tables, then we could do


~~~{.cxx}
    std::ofstream output("paper.input");
  	TIter next(allPlots);
  	TObject* o = 0;
  	Bool_t first = true;
  	while ((o = next())) {
	  Option_t* option = (first ? "hc" : "");
  	  if (o->IsA()->InheritsFrom(GraphSysErr::Class())) {
  	    GraphSysErr* g = static_cast<GraphSysErr*>(o);
  	    g->Export(output, option);
  	  }
  	  else if (o->IsA()->InheritsFrom(TSeqCollection::Class())) {
  	    TSeqCollection* c = static_cast<TSeqCollection*>(o);
  	    GraphSysErr::Export(c, output, option);
  	  }
      else
  	    // Ignore other stuff in collection.
  	  continue;
        first = false;
    }
		
~~~


### Import dataset from Durham database input file

If one has an _input_ formatted file of a Durham database dataset
(a single plot), one can import it into a GraphSysErr.

~~~{.cxx}
    GraphSysErr* g = GraphSysErr::Import("data.input", 1);
~~~

The second argument is the column number of the table in the file
_data.input_.  The systematic errors are called `sysX` where `X` is
the number (starting at 1).  To get the corresponding identifier one
can do

~~~{.cxx}
    UInt_t s1 = g->FindId(Form("sys%d",1));
    UInt_t s2 = g->FindId(Form("sys%d",2));
~~~
	
This ID can then be used to manipulate the systematic error.

To import multiple datasets, one do as done in the example
TestSample.

~~~{.cxx}
    TSeqCollection* c = GraphSysErr::Import("datasets.input");
~~~

Note, all tables are flattened.

## License

Copyright (c) 2014 Christian Holm Christensen <cholm@nbi.dk>

The code is licensed under the LGPL-3

<!-- Local Variables: -->
<!--   mode: markdown -->
<!-- End: -->
