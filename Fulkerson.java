import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Random;
public class Fulkerson{
    public static void fulkerson(List<ArrayList<Object>> graph,int source,int sink,List<Vertex> vertices,int algorithm_num,double[][] information,int graph_edges,int long_path_edge){
        List<ArrayList<Object>> residualgraph=new ArrayList<>(graph);
        double max_flow=0;
        int count_edges=0;
        int total_edges=0;
        int path=0;
        int num_path=0;
        while (true) {
            double cal_flow=0;
            
            
            switch (algorithm_num) {
                case 0:
                  //System.out.println("Running algorithm 1");
                  Object[] val=find_augmentingpath1(residualgraph,source,sink,vertices,count_edges,1);
                  cal_flow=(double)val[0];
                  count_edges=(int)val[1];
                 break;
                case 1:
                  //System.out.println("Running algorithm 2");
                  val=find_augmentingpath2(residualgraph,source,sink,vertices,count_edges,1);
                  cal_flow=(double)val[0];
                  count_edges=(int)val[1];
                  break;
                case 2:
                  //System.out.println("Running algorithm 3");
                  val=find_augmentingpath3(residualgraph,source,sink,vertices,count_edges,0);
                  cal_flow=(double)val[0];
                  count_edges=(int)val[1];
                  break;
                case 3:
                //System.out.println("Running algorithm 4");
                  val=find_augmentingpath4(residualgraph,source,sink,vertices,count_edges,1);
                  cal_flow=(double) val[0];
                  count_edges=(int) val[1];
                  break;
            }
            //double cal_flow=find_augmentingpath4(residualgraph,source,sink,vertices);
            if(cal_flow==0 || cal_flow==Double.POSITIVE_INFINITY || cal_flow==Double.NEGATIVE_INFINITY){
                break;
            }
            num_path=num_path+1;
            //System.out.println("Path flow: "+cal_flow);
            max_flow=max_flow+cal_flow;
            total_edges=total_edges+count_edges;
            //System.out.println("Max  flow: "+max_flow);
            
        }
        //System.out.println("Longest Path edge:"+long_path_edge);
        //System.out.println("Paths:"+num_path);
        //System.out.println("Total edges:"+total_edges);
        //System.out.println("Graph Edges:"+graph_edges);
        information[algorithm_num][0]=num_path;
        if(num_path==0){
            information[algorithm_num][1]=0;

        }
        else{
            information[algorithm_num][1]=total_edges/num_path;

        }
        
        information[algorithm_num][2] = Math.round((information[algorithm_num][1] / long_path_edge) * 1e6) / 1e6;
        information[algorithm_num][3]=graph_edges;
        information[algorithm_num][4]=max_flow;
        //System.out.println("Total flow"+max_flow);
        
    }
    public static Object[] find_augmentingpath1(List<ArrayList<Object>> residualgraph,int source,int sink,List<Vertex> vertices,int count_edges,int min_max_bit){
        double[] parent = new double[residualgraph.size()];
        Arrays.fill(parent, -1);
        double[] dist = new double[residualgraph.size()];
        boolean[] visited = new boolean[residualgraph.size()];
        int dest_idx;
        int flag=0;
        for (int i = 0; i < residualgraph.size(); i++) {
            dist[i] = Double.POSITIVE_INFINITY;
        }
        dist[source] = 0;
        double flow=0;

        PriorityQueue<Pair> q=new PriorityQueue<>(Comparator.comparingDouble(Pair -> Pair.path));
        q.add(new Pair(source, 0));
        parent[source]=source;
        while (!q.isEmpty()) {
            Pair curr=q.poll();
            // System.out.println("Currently exploring"+curr.n);
            // System.out.println("______");
            if(curr.n==sink){
                flag=1;
                break;
            }
            if (!visited[(int) curr.n]) {
                visited[(int) curr.n] = true;
        List<Edge> edges=(List<Edge>) residualgraph.get(curr.n).get(1);
        for(int i=0;i<edges.size();i++){
            Edge e=edges.get(i);
            int dest=vertices.indexOf(e.v);
        //    System.out.println("comparing Edges: "+curr.n+","+dest+","+e.cap);
            if(parent[dest]==-1 && e.cap>0){
                if(dist[curr.n]+e.cap<dist[dest]){
                   // System.out.println("Updating distance: "+e.cap);
                    dist[dest]=dist[curr.n]+e.cap;
            parent[dest]=curr.n;
            q.add(new Pair(dest,dist[dest]));
            }
        }
        }
        }
    }
            
      int num_edges=0;
            
        if(flag==1){
        Object[] val=add_flow(residualgraph,parent,source,sink,vertices,min_max_bit);
        flow=(double)val[0];
        num_edges=(int)val[1];
        if(flow>0){
            update_graph(residualgraph,parent,flow,source,sink,vertices);

        }

    }
        flag=0;
        return new Object[]{flow,num_edges};



    }
    public static Object[] find_augmentingpath2(List<ArrayList<Object>> residualgraph,int source,int sink,List<Vertex> vertices,int count_edges,int min_max_bit){
        double[] parent = new double[residualgraph.size()];
        int counter=1000;
        Arrays.fill(parent, -1);
        double[] dist = new double[residualgraph.size()];
        boolean[] visited = new boolean[residualgraph.size()];
        //int dest_idx;
        int flag=0;
        for (int i = 0; i < residualgraph.size(); i++) {
            dist[i] = Double.POSITIVE_INFINITY  ;
        }
        dist[source] = 0;
        double flow=0;

        PriorityQueue<Pair> q = new PriorityQueue<>(Comparator.comparingDouble(Pair -> Pair.path));
        q.add(new Pair(source,0));
        parent[source]=source;
        while (!q.isEmpty()) {
            Pair curr=q.remove();
            // System.out.println(curr.n);
            // System.out.println("______");
            if(curr.n==sink){
                flag=1;
                break;
            }
            if (!visited[(int) curr.n]) {
                visited[(int) curr.n] = true;
        List<Edge> edges=(List<Edge>) residualgraph.get(curr.n).get(1);
        for(int i=0;i<edges.size();i++){
            Edge e=edges.get(i);
            int dest=vertices.indexOf(e.v);
        //    System.out.println("comparing Edges: "+curr.n+","+dest+","+e.cap);
            
            if(e.cap>0){
                if(dist[dest]==Double.POSITIVE_INFINITY){
                    q.add(new Pair(dest, counter));
                    counter=counter-1;

            }
                
                if(dist[curr.n]+e.cap<dist[dest]){
                    // System.out.println("Selecting edge: "+e.cap);
                    
                    dist[dest]=dist[curr.n]+e.cap;
                    parent[dest]=curr.n;
            // System.out.println("Pushing values:"+dest+","+dist[dest]);

            q.add(new Pair(dest,dist[dest]));
            }
        }
        }
        }
    }
      int num_edges=0;
            
        if(flag==1){
        Object[] val=add_flow(residualgraph,parent,source,sink,vertices,min_max_bit);
        flow=(double)val[0];
        num_edges=(int)val[1];
        if(flow>0){
            update_graph(residualgraph,parent,flow,source,sink,vertices);

        }
    }
        flag=0;
        return new Object[]{flow,num_edges};



    }
    public static Object[] find_augmentingpath3(List<ArrayList<Object>> residualgraph,int source,int sink,List<Vertex> vertices,int count_edges,int min_max_bit){
        double[] parent = new double[residualgraph.size()];
        Arrays.fill(parent, -1);
        double[] dist = new double[residualgraph.size()];
        boolean[] visited = new boolean[residualgraph.size()];
        //int dest_idx;
        int flag=0;
        for (int i = 0; i < residualgraph.size(); i++) {
            dist[i] = Double.NEGATIVE_INFINITY;
        }
        dist[source] = 0;
        double flow=0;

        PriorityQueue<Pair> q = new PriorityQueue<>(Comparator.comparingDouble(Pair -> -Pair.path));
        q.add(new Pair(source, Double.POSITIVE_INFINITY));
        parent[source]=source;
        while (!q.isEmpty()) {
            Pair curr=q.remove();
            // System.out.println("sink is"+sink);
            // System.out.println(curr.n);
            // System.out.println("______");
            if(curr.n==sink){
                flag=1;
                break;
            }
            if (!visited[(int) curr.n]) {
                visited[(int) curr.n] = true;
        List<Edge> edges=(List<Edge>) residualgraph.get(curr.n).get(1);
        for(int i=0;i<edges.size();i++){
            Edge e=edges.get(i);
            int dest=vertices.indexOf(e.v);
            // System.out.println("comparing Edges: "+curr.n+","+dest+","+e.cap);
            if(parent[dest]==-1 && e.cap>0){
                
                if(dist[curr.n]+e.cap>dist[dest]){
                //  System.out.println("Selecting edge: "+e.cap);
                    
                    dist[dest]=dist[curr.n]+e.cap;
                    parent[dest]=curr.n;
                   
                    
                    q.add(new Pair(dest,dist[dest]));
            }
        }
        }
        }
    }
        int num_edges=0;    
    //   System.out.println("Flag is"+flag);
            
        if(flag==1){
        Object[] val=add_flow(residualgraph,parent,source,sink,vertices,min_max_bit);
        flow=(double)val[0];
        num_edges=(int)val[1];
        if(flow>0){
            update_graph(residualgraph,parent,flow,source,sink,vertices);

        }
    }
        flag=0;
        return new Object[]{flow,num_edges};



    }
    public static Object[] find_augmentingpath4(List<ArrayList<Object>> residualgraph,int source,int sink,List<Vertex> vertices,int count_edges,int min_max_bit){
        Random random = new Random();
        double[] parent = new double[residualgraph.size()];
        //int counter=1000;
        Arrays.fill(parent, -1);
        double[] dist = new double[residualgraph.size()];
        boolean[] visited = new boolean[residualgraph.size()];
        //int dest_idx;
        int flag=0;
        for (int i = 0; i < residualgraph.size(); i++) {
            dist[i] = Double.POSITIVE_INFINITY  ;
        }
        dist[source] = 0;
        double flow=0;

        PriorityQueue<Pair> q = new PriorityQueue<>();
        q.add(new Pair(source,0));
        parent[source]=source;
        while (!q.isEmpty()) {
            Pair curr=q.remove();
            // System.out.println(curr.n);
            // System.out.println("______");
            if(curr.n==sink){
                flag=1;
                break;
            }
            if (!visited[(int) curr.n]) {
                visited[(int) curr.n] = true;
        List<Edge> edges=(List<Edge>) residualgraph.get(curr.n).get(1);
        for(int i=0;i<edges.size();i++){
            Edge e=edges.get(i);
            int dest=vertices.indexOf(e.v);
            // System.out.println("comparing Edges: "+curr.n+","+dest+","+e.cap);
            
            if(parent[dest]==-1 && e.cap>0){
                if(dist[dest]==Double.POSITIVE_INFINITY){
                    int randomValue = random.nextInt(1000) + 1;
                    q.add(new Pair(dest, randomValue));
                    //counter=counter-1;

            }
                
                if(dist[curr.n]+e.cap<dist[dest]){
                // System.out.println("Selecting edge: "+e.cap);
                    
                    dist[dest]=dist[curr.n]+e.cap;
                    if(dest!=curr.n){
                        parent[dest]=curr.n;

                    }
                    // System.out.println("Pushing values:"+dest+","+dist[dest]);

            q.add(new Pair(dest,dist[dest]));
            
            }
        }
        }
        }
    }
       int num_edges=0;
            
        if(flag==1){
        Object[] val=add_flow(residualgraph,parent,source,sink,vertices,min_max_bit);
        flow=(double)val[0];
        num_edges=(int)val[1];
        if(flow>0){
            update_graph(residualgraph,parent,flow,source,sink,vertices);

        }
    }
        return new Object[]{flow,num_edges};



    }
    
    public static Object[] add_flow(List<ArrayList<Object>> residualgraph,double parent[],int source,int sink,List<Vertex> vertices,int min_max_bit){
        double minimum_cap=Double.POSITIVE_INFINITY;
        int num_edges=0;
        
        
        int curr=sink;
        
        while (curr != source) {
            int p = (int) parent[curr];
            
                // System.out.println("Current: " + curr);
                // System.out.println("Parent value: " + p);
                if (p == -1) {
                    System.out.println("Error: Source is unreachable!");
                    return new Object[]{-1, num_edges};
                }
                
                List<Edge> edges = (List<Edge>) residualgraph.get(p).get(1);
                for (int i = 0; i < edges.size(); i++) {
                    Edge e = edges.get(i);
                    if (e.v.equals(vertices.get(curr))) {
                        minimum_cap = Math.min(minimum_cap, e.cap);
                        num_edges = num_edges + 1;
                        break;
                    }
                }
                curr = p;
            } 
        
        

        return new Object[]{minimum_cap,num_edges};
}

    public static void update_graph(List<ArrayList<Object>> residualgraph,double parent[],double flow,int source,int sink,List<Vertex> vertices){
        int curr=sink;
        while (curr!=source) {
            int p=(int) parent[curr];
            List<Edge> edges=(List<Edge>)residualgraph.get(p).get(1);
            for(int i=0;i<edges.size();i++){
                Edge e=edges.get(i);
                if(e.v.equals(vertices.get(curr))){
                    e.cap=e.cap-flow;
                    break;
                }
            }
            
            
            curr=p;


    }
}
public static void helper_sap(List<ArrayList<Object>> graph,int source,int sink,List<Vertex> vertices,int algorithm_num,double[][] information,int graph_edges,int long_path_edge){
    List<ArrayList<Object>> cap1graph=new ArrayList<>(graph);
    for(int i=0;i<vertices.size();i++){
        // System.out.println("Vertex: "+cap1graph.get(i).get(0));
        for(int j=0;j<((List<Edge>) cap1graph.get(i).get(1)).size();j++){
            Edge e=((List<Edge>) cap1graph.get(i).get(1)).get(j);
            e.cap=1;
        }
    }
    fulkerson(cap1graph, source, sink, vertices,algorithm_num,information,graph_edges,long_path_edge);
    //dijkstra(cap1graph, source, sink, vertices);


    

}
}