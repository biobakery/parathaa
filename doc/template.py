#' **Feel free to edit/replace the template code below here and start working on the analysis.** 
#'  
#+ echo=False
# import pandas and matplotlib 
import pandas as pd 
import matplotlib.pyplot as plt 
import time
from anadama2 import PweaveDocument
document = PweaveDocument()

vars = document.get_vars()

#' Date: <%= time.strftime("%d, %b %Y") %>  
#'
#' # Demo Analysis Study Title
#'
#' # Introduction
#'  <% print(vars["introduction_text"]) %>Please follow the following code pattern while doing an analysis. 
#' This an example of a document that can be published using Pweave. Text is written using markdown (`#'`)
#' and python code between <%=%> are executed and results are included in the resulting pdf document.
#' You can define various options for code chunks to control code execution and formatting [see
#' Pweave docs](https://anadama2.readthedocs.io/en/latest/document.html).



#' # Examples
#' ### read_table example
#' 
#' ```
#' Description:
#' Read the table from a text file with the first line 
#' the column names and the first column the row names.
#'
#' Parameters:	
#' file (str) – The file to read
#' invert (bool) – Invert the table rows/columns after reading
#' delimiter (str) – The delimiter present in the file
#' only_data_columns (bool) – Remove the header and row names
#' format_data (function) – A function to use to format the data
#' ```
#' ##### Example Output:  
#+ echo=False
#'
#'
#' ### Displaying images from visualization modules-example
#' ![boxplots](../viz/boxplots.png)
#' The above boxplots is an example visualization output from plots.py
#' displayed as the markdown image in the Pweave pdf report. 
#'
#' ![boxplots](../viz/barplots.png)
#' The above bar plots is an example visualization output from plots.py
#' displayed as the markdown image in the Pweave pdf report.  

#'
#' # MISC Markdown Examples
#'
#' # h1 Heading 
#' ## h2 Heading
#' ### h3 Heading
#' #### h4 Heading
#' ##### h5 Heading
#' ###### h6 Heading

#' ## Tables

#' | Option | Description |
#' | ------ | ----------- |
#' | data   | path to data files to supply the data that will be passed into templates. |
#' | engine | engine to be used for processing templates. Handlebars is the default. |
#' | ext    | extension to be used for dest files. |

#' Right aligned columns

#' | Option | Description |
#' | ------:| -----------:|
#' | data   | path to data files to supply the data that will be passed into templates. |
#' | engine | engine to be used for processing templates. Handlebars is the default. |
#' | ext    | extension to be used for dest files. |


#' ## Emphasis

#' **This is bold text**
#' __This is bold text__
#' *This is italic text*
#' _This is italic text_
#' ~~Strikethrough~~