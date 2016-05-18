# Causal Structure Inference (CSI)

## The Purpose of the Algorithm

CSI is a network inference algorithm, designed to reconstruct regulatory models from time course data. The required input is a time course experiment, with the algorithm identifying relationships where one gene's expression profile explains another gene's expression profile at a later time. This accounts for all sorts of biological delays in the regulatory signal, such as translation and subsequently getting the regulating protein back into the nucleus. The delay between the trends in the regulator (parent) gene and the regulated (child) gene is set to be one time point of the time course. CSI was shown to outperform a number of other network inference algorithms when tested on synthetic DREAM4 data in [Penfold et al., 2011][penfold2011].

## Basic Input/Output

CSI accepts normally distributed data on input, so if you want to use RNA-Seq data on input consider log-transforming it beforehand. The output of CSI is a collection of CSV files with the edge probabilities in the network models, one obtained as the MAP estimate and one as the marginal estimate. These can be converted to GML files (which are compatible with a number of programs, including Cytoscape) with the hCSI_MarginalThreshold app. A webapp is also provided for interactive browsing of the resulting network model, which can be of use if the gene pool in the analysis was reasonably small (50 being an upper end estimate for reasonably clear visualisation).

## How Does It Work?

CSI tries to explain the expression of a child gene by taking the expression of potential parent genes at the previous time point. All possible parent combinations (from 0 up to a fixed maximum number of parents, due to computational constraints) are evaluated for any individual child gene. For each parent combination, a zero-mean Gaussian process is fit to the data, trying to explain the child's expression. The quality of any individual fit is measured in its likelihood, with Expectation-Maximisation being employed to help maximise the overall fit quality defined as the sum of the likelihoods of all the individual fits. Upon the completion of EM, the probability of a particular parent set being the "correct" one for the child gene is defined as the likelihood of that fit defined by the final likelihood sum. For more details, consult [Penfold et al., 2011][penfold2011].

## Test Run

If you want to take CSI out for a spin without using your own data, this can be done with the aid of one of the 10-gene synthetic networks originally used in the CSI and hCSI publications. The dataset to be used on input can be found at `ktpolanski/hcsi_testdata/dream4_5.csv` under Community Data. Leave all the parameter values as defaults.

## Input in Detail

**Note: CSI is quite computationally intensive.** Consider doing some preliminary analysis before feeding your data into it, with gene totals not exceeding the order of hundreds being preferable. Some good ways to reduce dimensionality include performing differential expression analysis and limiting the gene lists to transcription factors only.

### Gene Expression CSV

**Obligatory input.** Comma-delimited file, with expression data ordered to have genes as rows and time points as columns. In terms of headers, the first column should contain gene IDs, the first row should contain replicate names (repeated for each time point part of the replicate), and the second row should contain the corresponding time of the time point in that replicate. For reference on formatting, consult `ktpolanski/hcsi_testdata/dream4_5.csv` under Community Data. You can use multiple conditions on input, but CSI will treat them as replicates and not individual conditions, subsequently attempting to infer a joint regulatory model across all of the provided data.

### Parental Set Depth

**Default:** 2

When evaluating parental set combinations, a limitation is put on up to how many parents to sample from the parent pool to create the combinations. As this depth is increased, the number of parental sets to evaluate drastically goes up, making the problem less computationally tractable. Increasing this value above 2 is not recommended, especially if the dataset features hundreds of genes. At the same time, lowering this value to 1 will result in uninformative analysis results, missing out on a lot of combinatorial regulatory action.

### Gaussian Process Prior

**Default:** 10;0.1

Part of CSI is placing an assumption on how we expect the hyperparameter values to be distributed, and in [Penfold et al., 2011][penfold2011] the hyperparameters were assumed to be Gamma distributed with a shape parameter of 10 and a scale parameter of 0.1. If you wish to alter the shape/scale parameters of the gamma distribution, provide them as shape and scale with just a semicolon between them.

### Process Count

**Default:** 8

CSI features multiple levels of parallelisation, aiming to speed up the run time of the algorithm. This parameter controls how many threads are opened, and matches the default allocation of an iPlant node.

### Weight Truncation

**Default:** 1e-5

CSI's EM optimisation of the individual parent set fits sees a lot of computations of likelihoods, which is the computational bottleneck of the algorithm. If this option is provided, then some of the more suboptimal parent set combinations (whose probabilities constitute less than this fraction of the best parental set's probability) will be skipped to make the algorithm run quicker. Reducing this value will slightly improve result accuracy, at a computational cost. Increasing this value can lead to uninformative results as a high number of parental sets become dropped.

### No Data Standardisation

CSI fits are calibrated to normally distributed expression data with zero mean and unit variance. By default, CSI standardises the expression data it receives on a per-gene, per-condition basis to match this requirement. If you believe you have strong reason to not standardise the data automatically within the code, check this box.

### Weight Sampling

If checked, the starting values of the probabilities of the parental sets will be sampled from a Gamma distribution, with parental sets featuring fewer genes being preferred. If left unchecked, the starting values of the probabilities will be uniform across all of the parental sets.

## Output in Detail

### `csi-MAP.csv` / `csi-marginal.csv`

These are the regulatory network models, with parents as columns and children as rows, and each edge given a probability score (the higher the better). Introduce a threshold (for example, in the hCSI_MarginalThreshold app) to turn them into a standard, binary network. The MAP model is obtained by taking the single best parental set combination (the one with the highest probability) for each child, whilst the marginal model is obtained by summing the probabilities of all the parental sets featuring any individual parent. Typically, the marginal model is a bit more useful as it isn't limited to a single parental set fit in the information it provides.

### `csi_out.csv`

The raw CSI output, a comma-delimited file featuring the child gene in the first column, colon-delimited parent genes in the second column, the weight/probability of the parental set in the third column, and columns four to six being the hyperparameters of the optimised Gaussian process fits for that child gene. Provided in case more in depth information than that provided in the MAP/marginal models is desired.

### `csi_out.h5`

The same as the above, but provided in .h5 format. Used for preparing the webapp and computing the MAP/marginal models.

### `html/`

A webapp, allowing for an interactive visualisation and tuning of the nodes/edges present in the final exported network model. Open by double clicking on `csi.html` after downloading the results to your machine. Control the edge probability threshold in the Minimum Weight field, and switch between MAP/marginal models under the Network tab. The network visualisation itself is interactive and can have the nodes dragged around and the network zoomed in/out on. Select/deselect genes by clicking on them in the gene list, and see the expression profile of its chosen parents in the plot box below the space allocated for the network graph. Once content with a model, export it to a GML under the File tab. **Using the webapp is not recommended if more than 50 genes were used in the inference, as the model may become too large for the webapp to handle and visualise in an informative manner. GMLs can also be obtained through the hCSI_MarginalThreshold app.**

[penfold2011]: http://rsfs.royalsocietypublishing.org/content/1/6/857.short