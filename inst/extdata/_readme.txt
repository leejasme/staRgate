The files included in the extdata folder comprises of:
1. biexp_transf_parameters_x50.csv: 
    an example table of how to specify the parameters
    for biexponential transformation when pre-processing the GatingSet
    This is the suggested set of parameters as a starting point, the user
    may choose to modify as necessary based on the specific use case.
    a csv format must be supplied for the function getBiexpTransformGS() to work.

2. comp_mat_example_fcs.csv: 
    an example of compensation matrix. this corresponds to the 
    example_fcs.fcs also contained in this folder. 
    note that this is directly exported from flowJo as csv 
    a csv format must be supplied for the function getCompGS() to work.

3. example_fcs.fcs:
    this is an example FCS file included to run the tutorial vignette. 
    this FCS file was generated as a QC file at the Immune Monitoring Facility
    at Memorial Sloan Kettering Cancer Center on a BD FACSymphony.
    For illustration purposes, the origial FCS file was concatenated in flowJo
    to the first 30k events/cells acquired in order to reduce run time and file size

4. gating_template_x50_tcell.csv:
    this is a modified gating template based on the specified 29-marker panel. 
    the gating template is required to run the {openCyto} pre-gating step to 
    identify the key parent populations.
    for examples and more indepth explanations on how to modify the gating template, 
    please refer to the {openCyto} documentation at
    https://doi.org/doi:10.18129/B9.bioc.openCyto


5. pos_peak_thresholds.csv: 
    example and suggested values for positive peak thresholds for running the
    density gating. These parameters are obtained from a parameter tuning 
    exercise done on 1 dataset of 143 flow samples. 
    the user may choose to change the values based on the expectations for
    the specific use case, which may include different staining antibodies or
    different markers included on the flow panel.
