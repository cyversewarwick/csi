# Developer Notes for CSI #

I'll try to keep notes relevant to the translation of CSI to Python in
this document.

## Data from the DREAM challenge ##

Test data for CSI from the DREAM project is in `Input/Demo_DREAM.csv`.
It contains two header rows, the first defines the 5 "Treatments"
while the second row the 21 evenly spaced X values for the timeseries.
The first column contains the 10 gene names and the remaining 10 by
(21*5) cells contain the data values.

# From Chris #

## Parental Sets ##

Calculation of parental sets should probably be done outside of the
"main" CSI code—maybe in user friendly calling part.  Use case is to
optimise hyperparameters using a smaller set (i.e. truncate at two
parents) then run full algorithm with larger set (i.e. three parents).

# Displaying Results Interactively #

## Visual Elements ##

1. Graph/Network of genes as nodes and edges marginal likelihood over
   parental sets.  Hovering over nodes displays name, edge line weight
   proportional to marginals.

2. Plots of raw data; selected node black, genes having a casual
   effect on selected gene in one colour and genes our node has a
   causal effect on in another colour.  Line weight proportional to
   marginal likelihood.

3. A table of genes: name, number child and parent genes at current
   cutoff

4. Table of parental sets for a selected gene, highlighted by cutoff

5. Slider for cutoff on marginal likelihoods of parental sets.

6. Menubar for setting options, exporting data

## User Stories ##

I want to…

1. get a general idea of the network structure while easily cutting
   the graph at different marginals (uses: graph and cutoff)

2. make sure that CSI is "doing the right thing" by looking at
   expression profiles and comparing them to their parental sets (uses
   graph and plots)

3. look at my genes of interest and explore the network near them
   (uses graph and something else)

4. use Cytoscape to look at the network cut as appropriate.  I want
   the graph in a .SIF file and the edges in a .EDA file.

5. get citation information

# Files #

## Control Flow ##

`csi` is the entry point of the program, starts the GUI and responds
to button presses.  `panelfuncs` contains a list of panels to display,
all subdirectories are added to Matlab's search path.  The data is
loaded, parsed, genes and transcription factors selected and finally a
structure is passed to `run_CSI`.  This structure contains:

 * `filename` string
 * `data` a double matrix, as per demo data file
 * `genenames` a cell array of gene names
 * `genedesc` a cell array of gene descriptions
 * `replicatenames` a cell array of replicate names
 * `orig_startingCols` not sure
 * `startingCols` matrix describing which columns of the data refer to
   which treatment
 * `time_values` double vector, as per demo data
 * `toplot` scalar numeric value
 * `tf` double vector, indicies of genes to treat as transcription factors
 * `gene_idx` as per `tf`
 * `params`: has members
   *  `type` either `CSI` or `Hierarchical CSI`
   * `Pr` priors for hyper parameters
   * `inference` 1 for EM or 2 for MCMC (only for CSI).
   * `sparse` optimisation (only for EM)
   * `N` number of MC steps (only for MCMC or Hierarchical)
   * `temp` temp for Hierarchical
   * `fixtemp` boolean for Hierarchical?
   * `indegree` where to truncate the in-degree
   * `dirname` output directory
   * `parEnv` not sure, empty matrix!

`run_CSI` is still within the GUI, more investigation needed!  It also
documents the contents of `data`, which appear to be consistent with
what I inferred.  The first stage prints out a description of the
options chosen along with a summary of the data size.  Next the data,
`D`, is normalised (within each gene; mean=0, sd=1) and split into two
cell arrays `Xin` and `Yin` , each `nreps` long.  The data for gene
`ijk` and replicate `i` is extracted, 

`run_CSI` appears to have a big distinction between CSI and HCSI, file
handles are set up for preprocessing, running for a single gene, and
postprocessing.  Focusing on the _simple_ Vanilla EM CSI seems like a
good starting point, will come back to other (more complicated?)
versions later.

The preprocessing for CSI takes all the matrices from the nice cell
array and concatenates them back into a matrix!  It leaves a structure
with two elements `.X` and `.Y` defined.  There is a choice of whether
a serial or parallel run will be attempted, but both will loop through
every gene with the parallel version passing each gene off to a
separate task---could be interesting to measure how much time is spent
distributing data around the cluster, probably not much as the data
should be "small".

The meat of CSI is in `run_csi_for_one_gene` (in `run_CSI.m`), with a
big divide depending on whether sparse processing has been selected.
Starting with the non-sparse implementation, we calculate `KSTAR` then
the `Pa` (parental set?) and pass off to `CSI_EM_v2` (WTF is the
`_v2`, why isn't there just a "current" CSI EM algorithm
distributed---Chris says it's because of a lack of version control).

`CSI_EM_v2` is in `Functions/General`, this gets called from the GUI
with four parameters.  First half of the code seems to be decoding
parameters, probably for compatibility with different callers?
Parameters are then "randomly" initialised, and the parameters
optimised to minimise `CSIEngine_v2`.  The results are then summarised
and returned to the caller.

`CSIEngine_v2` is also in `Functions/General`, and seems to be about
as far as we need to go.  Various hacks to deal with varargs that
appear unnecessary, followed by lots of calls to the GPML library to
get marginal likelihoods.  There appear to be `_v3` versions of lots
of functions, after speaking with Chris, these are likely the versions
that came from Ahmed Shifaz and are only partially incorporated into
the code.
