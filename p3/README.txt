to compile:
g++ bioinformatics_p2.cpp -o UPGMA
to run:
./UPGMA
note: 
stdout is doing this weird wrap-around thing and I'm not sure why, you can still tell what 
my algorithm is doing and the newick output is correct but some of the output to tell you merging is happening
gets cut off before the next one starts and stdout replaced a few characters in the beginning of the
"Newick format: (...)" output with the trailing parentheses of the newick tree. I tried using printf instead
of cout and ran into the same issue, so I don't know why this is happening but it doesn't affect the 
correctness of the output.