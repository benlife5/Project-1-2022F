Constructing the graph:
first find all nodes in the graph by taking all unique (k-1)mers
for each kmer, create an edge from it's prefix to its suffic 7mer nodes

Constructing the genome:
find the Eulerian path in the debrujin graph 
    indegree and outdegree of a node must be equal for all nodes for this to work (balanced and strongly connected)
```
EulerianCycle(Graph)
    form a cycle Cycle by randomly walking in Graph (don't visit the same edge twice!)
    while there are unexplored edges in Graph
        select a node newStart in Cycle with still unexplored edges
        form Cycle’ by traversing Cycle (starting at newStart) and then randomly walking
        Cycle ← Cycle’
    return Cycle
```
In order to find a cycle from our debrujin graph, we need to balance the graph:
> a nearly balanced graph has an Eulerian path if and only if adding an edge between its unbalanced nodes makes the graph balanced and strongly connected.

Pick a node and randomly start walking it. Once stuck (a cycle is complete), pick a node in that cycle with unused edges and repeat.
    Stop once the cycle length == # of edges
