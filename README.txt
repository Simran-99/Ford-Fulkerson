GRAPH GENERATOR:
This Java program creates random graphs and conducts a comprehensive analysis of their properties using various algorithms for finding augmenting paths. The generated graphs are saved 
as text files for further examination. The algorithms implemented are part of the Ford-Fulkerson maximum flow algorithm.

HOW TO USE:
Compile the Java code and run the program using the following commands:

javac GraphGenerator.java
java GraphGenerator

The program generates graphs for diverse combinations of parameters, such as the number of vertices (n), radius (r), and upper capacity (upperCap). The resulting graphs are saved in 
separate text files, e.g., graph1.txt, graph2.txt, etc.

GENERATE GRAPH:

The generated graphs include information about vertices, edges, source, and sink. The graph files are formatted for easy readability and analysis.

ALGORITHM ANALYSIS:

Different algorithms are employed to find augmenting paths, and their properties are analyzed. The results are printed in the console, detailing the number of augmenting paths, average length, mean proportional length, total edges, and max flow for each algorithm.

Parameters
n: An array of integers representing the number of vertices in the graphs.
r: An array of doubles representing the radius of the graphs.
upperCap: An array of integers representing the upper capacity of edges in the graphs.
File Naming Convention
The generated graph files adhere to the naming convention: graph{counter}.txt, where {counter} denotes a sequential number.

RESULTS
The results of the algorithm analysis are displayed in the console, providing insights into the performance of each algorithm. The output includes columns for the number of augmenting paths, average length of augmenting paths, mean proportional length, total edges, and max flow for each algorithm.

ADDITIONAL INFORMATION
The program utilizes the Ford-Fulkerson algorithm with distinct augmenting path selection strategies.
