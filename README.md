## Adaptive Mixed Team Orienteering Problem with Time Windows

Iterated Local Search is consisted mainly of two components:  
1. Local search
1. Shake step

### Split Local Search
We will have a list that holds the unvisited nodes. At each revision, before we enter the Split Local Search (SLS) phase, 
we will assign each unvisited node to a Solution bucket based on some criteria.  

After we finish the SLS phase which constructs a Route plan for all Solution buckets, 
we will connnect the unvisited nodes again to be ready for the next revision.

### Shake Step
Shake procedure will be applied at a single solution with multiple walks.
