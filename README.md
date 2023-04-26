# Gibbs_sampling
A Python program to find motif sequences in a number of DNA sequemnces

Assume that the motif that we are searching for is 6 bp wide (W = 6)

1. Randomly remove one of the sequences 0-5
2. For the remaining 5 sequences, randomly choose a 6 bp substring (boxes above)
3. Calculate the score of this substring as a matrix (see next slide)
4. Calculate the match of this matrix to every position in the removed sequence
5. Sample the removed sequence in terms of the probability of the match at each position
6. Add the sampled site sequence to the list of site sequences
7. If the new site sequence list scores better, retain it
8. Randomly choose a new substring in the group of 5 sequences
9. Assign the identified substring in the removed sequence at the chosen substring, and add the removed sequence back with the other 5 sequences
10. Randomly choose another sequence to remove, and repeat from step 3
11. Continue for a set number of iterations of until convergence
12. Output the site sequence list

![image](https://user-images.githubusercontent.com/57168335/234695916-ed2213d6-1451-4def-9401-1ef6d0a9ed40.png)

