import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.ThreadLocalRandom;


public class Main {

    private static List<Node> nodesCovered = new ArrayList<>();
    private static List<Node> invitees = new ArrayList<>();
    private static List<Node> nodesByAddedValue = new ArrayList<>();

    private static HashSet noDupSet = new HashSet();
    private static List<Node> covered = new ArrayList<>();

    private static double roundThreeDecimals(double d) {
        DecimalFormat twoDForm = new DecimalFormat("#.###");
        return Double.valueOf(twoDForm.format(d));
    }

    private static void initializeRandomNetwork(Graph g, int numNodes, double p) throws NodeNotFoundException {
        for (int i = 0; i < numNodes; i++) {

            //int index, double p, double degree, double rewiringProb ,boolean invited, boolean attended, int uniqueNodesCovered
            g.addNode(new Node(i, Math.random(), 0, 0, false, false, 0));
        }
        for (Node s : g.getNodes()) {
            for (Node e : g.getNodes()) {
                if (!s.equals(e)) {
                    g.addEdge(p, s, e);
                }
            }
        }
    }

    private static void initializeSmallWorldNetwork(Graph g, int numNodes, double p, double rewiringProb) throws NodeNotFoundException {

        for (int i = 0; i < numNodes; i++) {
            double attendanceProb = Math.random();
            //double rewiringProb = Math.random();
            //int index, double p, double degree, boolean invited, boolean attended, int uniqueNodesCovered
            g.addNode(new Node(i, attendanceProb, 0.0, rewiringProb,
                    false, false, 0));
            //System.out.println("Node added: " + i + " , attendance p: " + roundThreeDecimals(attendanceProb)
            //        + " , rewiring p: " + roundThreeDecimals(rewiringProb));
        }

        for (Node s : g.getNodes()) {
            for (Node e : g.getNodes()) {
                if (!s.equals(e)) {
                    //System.out.println("----------------------------------");
                    g.addEdge(p, s, e);
                    //System.out.println("Edge added from: " + s.getIndex() + " to " + e.getIndex());
                }
            }
        }

        for (Node node : g.getNodes()) {
            for (Node neighbour : node.getNodesReachable()) {

                if (Math.random() < node.getRewiringProb()) {
                    int randomNode = ThreadLocalRandom.current().nextInt(0, g.getNodes().size());
                    Node m = g.getNodeByIndex(randomNode);
                    if (!node.equals(m)) {
                        g.removeEdge(node, neighbour);
                        g.addEdge(1.0, node, m);
                    }
                }
            }
        }
    }

    private static void scaleFree(Graph g, int numNodes, int initialNumNodes, double degree) throws NodeNotFoundException {
        for (int i = 0; i < initialNumNodes; i++) {
            g.addNode(new Node(i, Math.random(), degree, 0, false, false, 0));
        }
        for (Node s : g.getNodes()) {
            for (Node e : g.getNodes()) {
                if (!s.equals(e)) {
                    g.addEdge(degree, s, e);
                }
            }
        }

        for (int k = 0; k < numNodes - initialNumNodes; k++) {

            Node h = new Node(g.getNodes().size(), Math.random(), degree, 0, false, false, 0);
            g.addNode(h);

            for (Node node : g.getNodes()) {
                float d = ((float) node.getNeighboursSize()) / g.getEdges().size();
                node.setDegree(d);
            }

            Collections.sort(g.getNodes(), Comparator.comparingDouble(Node::getDegree));
            g.addEdge(1, g.getNodes().get(0), h);
        }

    }

    private static void initializeScaleFreeNetwork(Graph g, int numNodes, int initialNumNodes) throws NodeNotFoundException {

        /*
        The Algorithm:
            Input: Number of Nodes N;
                   Initial number of nodes m0;
                   Offset Exponent a;
                   Minimum degree 1 <= d <= m0.
            Output: scale-free multigraph G = ({0,....,N-1}, E).

            1) Add m0 nodes to G.
            2) Connect every node in G to every other node in G, i.e. create a complete graph.
            3) Create a new node i.
            4) Pick a node j uniformly at random from the graph G. Set P = (k(j)/k_tot)^a.
            5) Pick a real number R uniformly at random between 0 and 1.
            6) If P > R then add j to i's adjacency list.
            7) Repeat steps 4 - 6 until i has m nodes in its adjacency list.
            8) Add i to the adjacency list of each node in its adjacency list.
            9) Add i to to the graph.
            10) Repeat steps 3 - 9 until there are N nodes in the graph.
         */

        //1) Add m0 nodes to G.
        for (int i = 0; i < initialNumNodes; i++) {
            //int index, double p, double degree, double rewiringProb ,boolean invited, boolean attended, int uniqueNodesCovered
            g.addNode(new Node(i, Math.random(), 0, 0, false, false, 0));
        }
        //2) Connect every node in G to every other node in G, i.e. create a complete graph.
        for (Node s : g.getNodes()) {
            for (Node e : g.getNodes()) {
                if (!s.equals(e)) {
                    g.addEdge(1.0, s, e);
                }
            }
        }

        for (int k = 0; k < numNodes - initialNumNodes; k++) {

            //3) Create a new node i.

            //int index, double p, double degree, double rewiringProb ,boolean invited, boolean attended, int uniqueNodesCovered
            Node i = new Node(g.getNodes().size(), Math.random(), 0, 0, false, false, 0);
            //System.out.println("Node " + i.getIndex() + " added...");
            //4) Pick a node j uniformly at random from the graph G.
            // Set P = (k(j)/k_tot)^a.
            // k(j) = degree of j; number of nodes j is connected to
            // k_tot = total number of edges in G
            // a = offset exponent = 1.22
            for (int m = 0; m < initialNumNodes; m++) {

                int random = ThreadLocalRandom.current().nextInt(0, g.getNodes().size());
                Node j = g.getNodeByIndex(random);
                //System.out.println("------------------------------------");
                //System.out.println("Random node j: " + j.getIndex());
                //System.out.println("J's neighbour size: " + j.getNeighboursSize());
                //System.out.println("Nodes in G: " + g.getNodes().size());
                float degree = ((float) j.getNeighboursSize()) / g.getNodes().size();
                j.setDegree(degree);
                //System.out.println("Degree of j: " + roundThreeDecimals(j.getDegree()));
                //5) Pick a real number R uniformly at random between 0 and 1.
                double R = Math.random();
                //System.out.println("Random number: " + R);
                //System.out.println("------------------------------------");
                //6) If P > R then add j to i's adjacency list.
                if (j.getDegree() >= R) {
                    //System.out.println("Edge added from " + i.getIndex() + " to: " + j.getIndex());
                    if (i.getEdgeTo(j) == null) i.addEdge(new Edge(i, j));
                }
            }
            //8) Add i to the adjacency list of each node in its adjacency list.
            for (Node node : i.getNodesReachable()) {
                node.addEdge(new Edge(node, i));
            }
            //9) Add i to to the graph.
            g.addNode(i);
        }
    }

    private static void randomHeuristic(int numInvitees, Graph g) {
        for (int i = 0; i < numInvitees; i++) {
            int random = ThreadLocalRandom.current().nextInt(0, g.getNodes().size());
            Node node = g.getNodeByIndex(random);
            if (!node.isInvited()) {
                node.setInvited(true);
            } else {
                numInvitees += 1;
                //System.out.println("Iterations changed: " + numInvitees);
            }
        }
    }

    private static void order(List<Node> nodes) {
        Collections.sort(nodes, Comparator.comparingDouble(Node::distribution));
        Collections.reverse(nodes);
    }

    private static void adaptiveGreedyStatic(int numInvitees, List<Node> nodes) {
        Node first = nodes.get(0);
        first.setInvited(true);
        nodesCovered.addAll(first.getNodesReachable());
        maximize(nodes);

        int counter = 0;
        while (counter < numInvitees - 1) {
            Node n = nodesByAddedValue.get(counter);
            if (!n.isInvited()) {
                n.setInvited(true);
                nodesCovered.addAll(first.getNodesReachable());
                maximize(nodes);
                counter++;
            } else {
                numInvitees++;
                counter++;
            }
        }
    }

    private static void maximize(List<Node> nodes) {
        nodesByAddedValue.clear();
        for (Node node : nodes) {
            if (!node.isInvited()) {
                int intersection = intersection(nodesCovered, node.getNodesReachable());
                node.setUniqueNodesCovered(node.getNeighboursSize() - intersection);
                nodesByAddedValue.add(node);
            }
        }
        Collections.sort(nodesByAddedValue, Comparator.comparingDouble(Node::getAddedValue).reversed());

    }

    //TODO: algorithm is: sort -> mark as invited -> if P > attendance prob -> marks as covered ->
    //TODO: -> add neighbours to covered set -> sort again to get next best nodes with biggest (neighbour size * attendance prob)
    //TODO: repeat for N rounds

    private static void adaptiveGreedyDynamic(int numInvitees, List<Node> nodes) {
        Node first = nodes.get(0);
        first.setInvited(true);
        if (first.getP() >= Math.random()) {
            first.setAttended(true);
            for (Node node : first.getNodesReachable()) {
                node.setCovered(true);
                nodesCovered.add(node);
            }
            first.setUniqueNodesCovered(first.getNeighboursSize());
        }
        maximize(nodes);

        int counter = 0;
        while (counter < numInvitees - 1) {
            Node n = nodesByAddedValue.get(counter);
            if (!n.isInvited()) {
                n.setInvited(true);
                if (n.getP() >= Math.random()) {
                    n.setAttended(true);
                    for (Node node : n.getNodesReachable()) {
                        node.setCovered(true);
                        nodesCovered.add(node);
                    }
                    maximize(nodes);
                }
                counter++;
            } else {
                numInvitees++;
                counter++;
            }
        }
    }

    private static void bestNextNodeHeuristic(int numInvitees, List<Node> nodes) {
        int counter = 0;
        Node bestFirst = nodes.get(0);
        //System.out.println("First node: " + bestFirst.getIndex() + ", prob: " + roundThreeDecimals(bestFirst.getP()) + ", neighbours: " + bestFirst.getNodesConnected());
        bestFirst.setInvited(true);
        invitees.add(bestFirst);
        while (counter < numInvitees-1) {
            for (Node node : nodes) {
                if (!invitees.contains(node)) {

                    int coverageForNodeInList = 0;
                    double averageCoverageForNodeInList = 0;
                    int sizeList = 0;
                    int differenceList = 0;
                    for (int i = 0; i < 100; i++) {

                        for (Node invitee : invitees) {
                            if (invitee.getP() >= Math.random()) {
                                nodesCovered.addAll(invitee.getNodesReachable());
                            }
                        }

                        sizeList = nodesCovered.size();

                        if (node.getP() >= Math.random()) {
                            nodesCovered.addAll(node.getNodesReachable());
                        }

                        differenceList = nodesCovered.size() - sizeList;
                    }


                    //System.out.println("Nodes covered before node " + node.getIndex() + " : " + sizeList + ", nodes covered size after: " + nodesCovered.size());
                    coverageForNodeInList += differenceList;
                    //System.out.println("Coverage for node in list: " + differenceList);

                    averageCoverageForNodeInList = (((double) coverageForNodeInList) / 100);
                    //System.out.println("Average coverage for node: " + node.getIndex() + ": " + averageCoverageForNodeInList);
                    node.setAverageCoverage(averageCoverageForNodeInList);
                    //System.out.println("Best coverage for node: " + node.getBestAverageCoverage());
                    //System.out.println("---------------------------------------------");
                }
                nodesCovered.clear();
            }
            Collections.sort(nodes, Comparator.comparingDouble(Node::getBestAverageCoverage).reversed());
            //for (Node node : nodes) {
            //    System.out.println("Node " + node.getIndex() + ", prob: " + node.getP() + " , neighbours: "
            //            + node.getNeighboursSize() + " , average coverage: " + roundThreeDecimals(node.getBestAverageCoverage()));
            //}
            for (Node n : nodes) {
                if (!invitees.contains(n)) {
                    n.setInvited(true);
                    invitees.add(n);
                    break;
                }
            }
            counter++;
            nodesCovered.clear();
        }
    }

    private static void greedyHeuristic(int numInvitees, List<Node> nodes) {
        int counter = 0;
        while (counter < numInvitees) {
            Node node = nodes.get(counter);
            if (!(node.isInvited()) && !(node.isCovered())) {
                node.setInvited(true);
                counter++;
            } else {
                numInvitees += 1;
                //System.out.println("Iterations changed: " + numInvitees);
                counter++;
            }
        }
    }

    private static void runTRounds(int rounds, int numInvitees, Graph g) {
        for (int i = 0; i < rounds; i++) {
            System.out.println("ROUND " + i);
            adaptiveGreedyStatic(numInvitees, g.getNodes());
        }
    }

    public static Integer intersection(List<Node> list1, List<Node> list2) {
        List<Node> list = new ArrayList<>();
        for (Node n : list1) {
            if (list2.contains(n)) {
                list.add(n);
            }
        }
        return list.size();
    }

    private static PrintWriter createFileI(int rounds, String typeGraph, String typeHeuristic, int numNodes, int initialNumNodes, int numInvitees) {
        PrintWriter pw = null;
        try {
            pw = new PrintWriter(new File(
                    "rounds=" + rounds + "-graph=" + typeGraph + "-heuristic=" + typeHeuristic + "-nodes=" +
                            numNodes + "-init=" + initialNumNodes + "-invited=" + numInvitees + ".csv"));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        return pw;
    }


    private static PrintWriter createFileF(int rounds, String typeGraph, String typeHeuristic, int numNodes, double probability, int numInvitees) {
        PrintWriter pw = null;
        try {
            pw = new PrintWriter(new File(
                    "rounds=" + rounds + "-graph=" + typeGraph + "-heuristic=" + typeHeuristic + "-nodes=" +
                            numNodes + "-prob=" + probability + "-invited=" + numInvitees + ".csv"));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        return pw;
    }

    public static void main(String[] args) {

        int averageAttended = 0;
        int averageUnique = 0;
        int averageTotal = 0;

        int numNodes = 200;
        double prob = 0.5;
        int numInv = 10;

        PrintWriter pw = createFileF(100, "smallWorld", "bestNext", numNodes, prob, numInv);

        long startTime = System.nanoTime();

        for (int i = 0; i < 100; i++) {
//
            Graph graph = new Graph();

            try {
                initializeSmallWorldNetwork(graph, numNodes, 0.1 ,prob);
            } catch (NodeNotFoundException e) {
                e.printStackTrace();
            }
            order(graph.getNodes());
//
            bestNextNodeHeuristic(numInv, graph.getNodes());
//
            //System.out.println("-------------------------------------------------");
//
//
            int countInvited = 0;
            int countAttended = 0;
            for (Node node : graph.getNodes()) {
                if (node.isInvited()) {
                    countInvited++;
                    if (node.getP() >= Math.random()) {
                        countAttended++;
                        covered.addAll(node.getNodesReachable());
                        noDupSet.addAll(node.getNodesReachable());
                    }
                }
            }
            //System.out.println();
//
            //System.out.println("Nodes invited:           " + countInvited);
            //System.out.println("Nodes attended:          " + countAttended);
            //System.out.println("Nodes covered in total:  " + covered.size());
            //System.out.println("Unique nodes covered:    " + noDupSet.size());
//
            //System.out.println("--------------------------------------------------");
//
            StringBuilder sb = new StringBuilder();
            sb.append("Round:");
            sb.append(',');
            sb.append(i);
            sb.append('\n');
            sb.append("Attended");
            sb.append(',');
            sb.append(countAttended);
            sb.append('\n');
            sb.append("Total covered");
            sb.append(',');
            sb.append(covered.size());
            sb.append('\n');
            sb.append("Unique covered");
            sb.append(',');
            sb.append(noDupSet.size());
            sb.append('\n');
            sb.append('\n');
            pw.write(sb.toString());

            averageAttended += countAttended;
            averageUnique += noDupSet.size();
            averageTotal += covered.size();

            covered.clear();
            noDupSet.clear();
        }

        pw.close();

        System.out.println("Average attended: " + ((float) averageAttended / 100));
        System.out.println("Average unique:   " + ((float) averageUnique / 100));
        System.out.println("Average total:    " + ((float) averageTotal / 100));

        long endTime = System.nanoTime();
        long totalTime = endTime - startTime;
        System.out.println("Time:             " + totalTime / 1000000000 + " secs");
    }

}