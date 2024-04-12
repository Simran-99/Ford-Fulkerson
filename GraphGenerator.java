import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Random;
import java.util.PriorityQueue;
//import path.to.Fulkerson;

class Vertex {
    double x, y;
    int n;

    Vertex(int n,double x, double y) {
        this.n=n;
        this.x = x;
        this.y = y;
    }
}

class Edge {
    Vertex u, v;
    double cap;

    Edge(Vertex u, Vertex v, double cap) {
        this.u = u;
        this.v = v;
        this.cap = cap;
    }
}

class Pair implements Comparable<Pair> {
    int n;
    double path;

    Pair(int n, double path) {
        this.n = n;
        this.path = path;
    }

    @Override
    public int compareTo(Pair p2) {
        return Double.compare(this.path, p2.path);
    }

}
public class GraphGenerator {   
    public static void main(String[] args) {
        int n_val[]={100,200};
        double r_val[]={0.2,0.3,0.15};
        int upper_cap[]={2,50,10};
        int counter=1;
        
        for(int i=0;i<n_val.length;i++){
            for(int j=0;j<r_val.length;j++){
                for(int k=0;k<upper_cap.length;k++){
                    String fname="graph"+counter+".txt";
                    System.out.println("n: "+n_val[i]+"r: "+r_val[j]+"cap: "+upper_cap[k]);
                    main_code(n_val[i],r_val[j], upper_cap[k],fname);
                    counter=counter+1;

                }
            }
        }
        
        

    }
    public static void main_code(int n, double r,int upperCap,String fname){
        // int n =n;
        // double r = 0.2;
        // int upperCap = 50;
        List<Vertex> vertices = generateVertices(n);
        // for (int i = 0; i < vertices.size(); i++) {
        //     Vertex vertex = vertices.get(i);
        //     //System.out.println("("+ i+","+ vertex.x + "," + vertex.y + ")");
        // }
        List<Edge> edges = generateEdges(vertices, r, upperCap);
        // for (int i = 0; i < edges.size(); i++) {
        //     Edge e = edges.get(i);
        //    // System.out.println("("+e.u.n+"," + e.u.x + "," + e.u.y + ") - ("+e.v.n+","   + e.v.x + "," + e.v.y + ") Capacity: " + e.cap);
        // }
        List<ArrayList<Object>> graph = generateGraph(vertices, edges);
        findSourceSink(vertices, edges, graph,edges.size(),fname);
    }
    private static List<Vertex> generateVertices(int n) {
        List<Vertex> vertices = new ArrayList<>();
        Random rand = new Random();

        for (int i = 0; i < n; i++) {
            double x = rand.nextDouble();
            double y = rand.nextDouble();
            vertices.add(new Vertex(i,x, y));
        }

        return vertices;
    }

    private static List<Edge> generateEdges(List<Vertex> vertices, double r, int upperCap) {
        List<Edge> edges = new ArrayList<>();
        Random rand = new Random();
        for (int i = 0; i < vertices.size(); i++) {
            for (int j = 0; j < vertices.size(); j++) {
                Vertex u = vertices.get(i);
                Vertex v = vertices.get(j);
                if(u!=v){
                if (u != v && (((u.x - v.x) * (u.x - v.x)) + ((u.y - v.y) * (u.y - v.y)) <= (r * r))) {
                    double val = rand.nextDouble();
                    if ((val < 0.5) && (!hasEdge(edges, u, v))) {
                        
                            edges.add(new Edge(u, v, rand.nextInt(upperCap - 1) + 1));
                        }
                    }
                    else {
                            if (!hasEdge(edges, v,u)) {
                                edges.add(new Edge(v, u, rand.nextInt(upperCap - 1) + 1));
                            }
                        }
                }
            }
            }
            return edges;
        }
        
    

    private static boolean hasEdge(List<Edge> edges, Vertex u, Vertex v) {
        for (Edge edge : edges) {
            if ((edge.u == u && edge.v == v) || (edge.u == v && edge.v == u)) {
                return true;
            }
        }
        return false;
    }

    private static List<ArrayList<Object>> generateGraph(List<Vertex> vertices, List<Edge> edges) {
        List<ArrayList<Object>> graph = new ArrayList<>();
        for (int i = 0; i < vertices.size(); i++) {
            ArrayList<Object> v = new ArrayList<>();
            v.add(vertices.get(i));
            List<Edge> e = new ArrayList<>();
            for (int j = 0; j < edges.size(); j++) {
                if (edges.get(j).u == vertices.get(i)) {
                    e.add(edges.get(j));
                }
            }
            v.add(e);
            graph.add(v);
        }
        return graph;
    }
    private static Object[] writegraph(List<Vertex> vertices, List<Edge> edges, List<ArrayList<Object>> graph,Vertex source,Vertex sink,String fname) {
   

        try (FileWriter fw = new FileWriter(fname)) {
            
            for (ArrayList<Object> vertexInfo : graph) {
                Vertex vertex = (Vertex) vertexInfo.get(0);
                fw.write("Vertex: (" + vertex.n + "," + vertex.x + "," + vertex.y + ")\n");
    
                List<Edge> edgesList = (List<Edge>) vertexInfo.get(1);
                for (Edge edge : edgesList) {
                    fw.write("Edge: (" + edge.u.n + "," + edge.u.x + "," + edge.u.y +
                             ") - (" + edge.v.n + "," + edge.v.x + "," + edge.v.y +
                             ") Capacity: " + edge.cap + "\n");
                }
            }
    
            // Write source and sink information to the file
            fw.write("Source: (" + source.n + "," + source.x + "," + source.y + ")\n");
            fw.write("Sink: (" + sink.n + "," + sink.x + "," + sink.y + ")\n");
        } catch (IOException e) {
            e.printStackTrace();
        }
        return new Object[]{source, sink, graph};
    
       
    }
    private static void findSourceSink(List<Vertex> vertices, List<Edge> edges, List<ArrayList<Object>> graph,int edges_size,String fname) {
    double information[][]=new double[4][5];
    Queue<Integer> q = new LinkedList<>();
    boolean[] vis = new boolean[vertices.size()];
    Random rand = new Random();
    int distance[]=new int[vertices.size()];
    int ridx = rand.nextInt(vertices.size());
    Vertex source = vertices.get(ridx);
    int long_path_edge=0;
    //System.out.println("Source: ("+source.n+","+ source.x + "," + source.y + ")");
    
    q.add(ridx);
    while (!q.isEmpty()) {
        int currIdx = q.remove();
        if (!vis[currIdx]) {
            vis[currIdx] = true;
        }
        
        List<Edge> edgesForCurr = (List<Edge>) graph.get(currIdx).get(1);
        for (int i = 0; i < edgesForCurr.size(); i++) {
            Edge e = edgesForCurr.get(i);
            if (!vis[vertices.indexOf(e.v)]) {
                q.add(vertices.indexOf(e.v));
                
            }
            distance[vertices.indexOf(e.v)] = distance[currIdx] + 1;
            List<Edge> edgesForCurrAdj = (List<Edge>) graph.get(vertices.indexOf(e.v)).get(1);

            for (int j = 0; j < edgesForCurrAdj.size(); j++) {
            Edge e1 = edgesForCurrAdj.get(j);
            distance[vertices.indexOf(e1.v)] = distance[vertices.indexOf(e.v)] + 1;
            
            
        }
    }
}
        int sink_idx=find_max(distance);
        Vertex sink = vertices.get(sink_idx);
    //     for(int i=0;i<distance.length;i++){
    //     System.out.println(distance[i]);
    // }
        long_path_edge=distance[sink_idx];
        
        
        if(source!=sink){
            //System.out.println("Sink: ("+sink.n+", "+ sink.x + "," + sink.y + ")");
            writegraph(vertices,edges,graph,source,sink,fname);
            Object[] sourceSinkInfo =writegraph(vertices,edges,graph,source,sink,fname);
                Vertex source1 = (Vertex) sourceSinkInfo[0];
                Vertex sink1 = (Vertex) sourceSinkInfo[1];
                List<ArrayList<Object>> graphInfo = (List<ArrayList<Object>>) sourceSinkInfo[2];
        
                for(int i=0;i<4;i++){
                    
                    List<ArrayList<Object>> graphz=deepCopy(graphInfo);
        
                    if(i==0){
                        Fulkerson.helper_sap(graphz, source1.n, sink1.n, vertices, i,information,edges_size,long_path_edge);
                    }
                    else{
                    Fulkerson.fulkerson(graphz, source1.n, sink1.n, vertices,i,information,edges_size,long_path_edge);
                    }
        
            }
            System.out.println("RESULTS:");
            System.out.println("Column 1: Number of Augmenting paths");
            System.out.println("Column 2: Average length of augmenting paths");
            System.out.println("Column 3: Mean proportional length");
            System.out.println("Column 4: Total edges");
            System.out.println("Column 5: Max Flow");
                for(int i=0;i<4;i++){
                for(int j=0;j<5;j++){
                    System.out.print(information[i][j]+" ");
                    
                }
                System.out.println();
            }

            }
            else{
                findSourceSink(vertices, edges, graph, edges_size,fname);
            }
            
    
    }
    private static List<ArrayList<Object>> deepCopy(List<ArrayList<Object>> graph){
        List<ArrayList<Object>> copied_graph=new ArrayList<>();
        for (ArrayList<Object> nodeinfo : graph) {
            Vertex v=(Vertex) nodeinfo.get(0);
            List<Edge> edgelist=(List<Edge>) nodeinfo.get(1);
            ArrayList<Object> n=new ArrayList<>();
            n.add(new Vertex(v.n, v.x, v.y));
            List<Edge> new_edge_list = new ArrayList<>();
            for(int i=0;i<edgelist.size();i++){
                Edge e=edgelist.get(i);
                new_edge_list.add(new Edge(e.u, e.v, e.cap));
            }
            n.add(new_edge_list);
            copied_graph.add(n);
        }
        return copied_graph;
    
    
    }
    private static int find_max(int arr[]) {
        int max_idx = 0;
    
        for (int i = 0; i < arr.length; i++) {
            if (arr[i] > arr[max_idx]) {
                max_idx = i;
            }
        }
    
        return max_idx;
    }
   
}