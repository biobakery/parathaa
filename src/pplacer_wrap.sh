#!/bin/bash

# references = $1
# queries = $2
# threads = $3


#99026
nseqs=$(wc -l $2 | cut -f 1 -d ".")

echo $nseqs
#determine number of chunks to divide into

((result = ($nseqs + $3 - 1)/ $3))

csplit $result

