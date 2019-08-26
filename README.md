# Target-RNASEQ-Analysis-GUI
A bokeh-based Graphical User Interface which is intended as a companion analysis tool for the Target RNA SEQ Protocol made by Grace Johnson and Darren Parker.

For now, in order to use this GUI you will have to download this folder and then navigate to this folder in command line and call:

``` bokeh serve --show folder_name ```

In the near future I want to push this to pip so that it can just be downloaded via pip.

## Tabs in the GUI

#### Control Panel Tab

When you open the GUI this is the first tab which opens. The goal of this software is to navigate you through the entire data analysis pipeline, started from the raw fastq files from the sequencer.  It will also save intermediary files, so that if you have done this analysis before you can skip to plotting the samples.  Here is what this tab looks like:

![alt text][cp_tab] 

Each of the columns of widgets performs a seperate function.  The first column automates the data-prepocessing needed before the reads can be aligned.  It removes all poly-A tails from the reads (which would cause them to mis-align) and also sorts all reads based on the barcodes that you provide.  Running this will automatically save the sorted and stripped files, and put the destination to this folder as the input to the next column. 

The next column handles the alignment of these reads back to the genome. All it needs is the file path to where the sorted and stripped files are kept.  Once again the aligned data file will automatically be written and saved. 

The last functional column will plot the aligned reads.  It requires the file output from the previous column in addition to the information about the primers used. Once you click plot here, the data for the other two tabs will be updated to reflect the change.

The last column on the far left is intended to be a help column.  There are a few FAQ that are included in the drop down.  Select the one you are interested in and press the *Help?* button to get the response. 

#### Primer Comparison Tab

This tab contains a log log plot that allows you to compare the relative prevalence of every type of primer across all of your experiments (experimental identity is determined by the barcodes of the reads).  To change which experiments you are comparing select a different subset from the dropdown menus.  In addition, you can select a specific subset of primers to highlight and label from the multi-select widget with all the primer names.  After making either one of these changes just select the *Plot* button and the changes should display. *In the near future I plan to add an option to control and change the lines which are displayed, but currently it just allows you to display a line at 1x, 2x, and 0.5x concentrations.* 

![alt text][pc_tab] 

#### Length Distribution Tab

![alt text][ld_tab] 

[cp_tab]: Images/control_panel_tab.PNG "Control Panel/Landing Page"
[pc_tab]: Images/prim_comp_tab.PNG "Comparison of all Primers Tab"
[ld_tab]: Images/len_comp_tab.PNG "Length Distribution for Given Primer"
