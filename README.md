# lakeMetabolism
Code and sample data to estimate lake ecosystem metabolism (pelagic primary production and respiration) from dissolved oxygen and other data. From Solomon et al. 2013 Limnology and Oceanography, DOI:10.4319/lo.2013.58.3.0849.

Hi there-

Some notes:

You can run everything by running the ‘script Annie.R’ file. This calls the main workhorse function ‘metabFunc_v6.r’, the function that really does the fitting ‘metabLoss_v4.r’, and other stuff.

There is a lot of extra code in here that was aimed entirely at the problem of massaging data sets from many different lakes into a common format. You may be able to do away with a lot of it.

The sample data that is in here is not ‘mine’, but comes from Lake Annie and the Archbold Biological Station. It is ok to play around with it to understand the code but please don’t use it for anything else.

A lot of time and effort went into coding up this model. We are happy to share it but if you use it we would appreciate being cited. The relevant paper is Solomon et al. 2013 Limnology and Oceanography.

The input data files have to be in a particular format (including column names etc) – look at the Annie data to see this

Thanks, and enjoy!
Chris
