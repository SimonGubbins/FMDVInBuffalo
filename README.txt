This folder contains the Matlab scripts and functions, as well the necessary data files, to implement the parameter
estimation methods in:

Jolles et al. (2020) Endemic persistence of a highly contagious pathogen: foot-and-mouth disease in its wildlife
reservoir. Science (in press)

DOI for code: 10.5281/zenodo.5121203

MATLAB REQUIREMENTS
The scripts/functions were run using Matlab version 2019b and require the Statistics and Machine Learning and
Parallel Computing toolboxes. However, they can be easily adapted to run without the Parallel Computing toolbox
by changing the "parfor" loop in the ParEst_[parameter] functions to a "for" loop.

DURATION OF MATERNAL ANTIBODIES (Supplementary Information S1)
ParEst_MaternalAbDuration.m - loads the data, implements the adaptive Metropolis scheme for each model/parameterisation
                              and computes the DIC
Lhood_MaternalAbDuration.m - computes the log likelihood and prior for the input parameters

MaternalAbDuration_Age0.txt - age of last protective titre for each calf (row) for each serotype (column)
MaternalAbDuration_Age1.txt - age of first non-protective titre for each calf (row) for each serotype (column)

ACUTE TRANSMISSION (Supplementary Information S2)
ParEst_AcuteTransmission.m - loads the data and implements the adaptive Metropolis scheme for each
                             model/parameterisation
Lhood_AcuteTransmission.m - computes the log likelihood and prior for the input parameters

LOOModelComparison_AcuteTransmission.m - implements the PSIS-LOO model comparison
PointwiseLhood_AcuteTransmission.m - computes the contribution to log likelihood for the input parameters for each
                                     buffalo (i.e. without combining them into a single log likelihood)

BuffaloTransmissionData.txt - transmission data; columns are:
                              1 - serotype (1-3 for SAT1-3)
                              2 - infection type (1-inoculated; 2-contact)
                              3 - time of last negative PCR result
                              4 - time of first positive PCR result
                              5 - time of last positive PCR result
                              6 - time of first subsequent negative PCR result

CARRIER TRANSMISSION (Supplementary Information S3)
ParEst_CarrierTransmission.m - contains the data, implements the adaptive Metropolis scheme for each
                               model/parameterisation and computes the DIC
Lhood_CarrierTransmission.m - computes the log likelihood and prior for the input parameters

CARRIER DURATION (Supplementary Information S4)
ParEst_CarrierDuration.m - contains the data, implements the adaptive Metropolis scheme for each
                           model/parameterisation and computes the DIC
Lhood_CarrierDuration.m - computes the log likelihood and prior for the input parameters


