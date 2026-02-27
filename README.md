# mhtexp2r
An attempt to port mhtexp2 from Stata to R.

The port was assisted greatly by Claude Opus 4.5, Opus 4.6, Sonnet 4.6, ChatGPT 5.2 and ChatGPT 5.

The code is provided as is without any guarantees. Use it, but you are responsible for what it does or does not do.

# mhtexp2 - Original README --The email address is not mine, but that of the original author
Implements the procedure detailed in List, Shaikh, and Vayalinkal (2023) (pdf: http://individual.utoronto.ca/atom/files/mht2.pdf). 

Please contact me (atom.vayalinkal@mail.utoronto.ca) with any questions/comments - thanks!


<h1> Installation instructions for Stata</h1>

This package extends the Stata `mhtexp` package, and can be installed by placing the mhtexp2.ado and mhtexp2.sthlp files in your stata installation’s "personal ado" folder. The location of your personal ado folder can be found using the stata command `personal`; for additional information, see here: https://www.stata.com/manuals13/u17.pdf#u17.7

Once installed, the package can be used by running the following command, from any do-file or console:
 
`mhtexp2 y, treatment(treatment_var) controls(control_var1 control_var2)`
 
Other than the new "controls" option, the options available are generally the same as the original mhtexp (available on ssc as `mhtexp`, and also at https://github.com/seidelj/mht). You also now have the option to use the procedure with and without studentization, and to pregenerate bootstrap samples (useful when running evaluating several different families using the same dataset). Please see the replication files for an example of how these operations work.
