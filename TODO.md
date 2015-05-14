Immediate

- Add content to the vignette stub

Long term

- Separate classes for 'DataFrame' and '\*Ranges' rowData

Possibilities?

- SummarizedExperiment virtual base class with derived classes

  - SummarizedExperimentDF
    @rowData: DataFrame

  - SummarizedExperimentGR
    @rowRanges: *Ranges; rowData() == mcols(rowRanges())

- SummarizedExperiment as 'DataFrame' base class,
  SummarizedExperimentGR as derived class

  - SummarizedExperiment
    @rowData: DataFrame

  - SummarizedExperimentGR
    @rowRanges: \*Ranges; no mcols() (rowData() from inheritted
    @rowData)

