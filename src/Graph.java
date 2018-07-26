import java.util.LinkedList;
import java.util.List;

public class Graph {

    private List<Node> nodes;
    private List<Edge> edges;

    public Graph() {
        this.nodes = new LinkedList<>();
        this.edges = new LinkedList<>();
    }

    public boolean addNode(Node node) {
        boolean isAdded = false;
        if (!nodes.contains(node)) {
            isAdded = this.nodes.add(node);
        }
        return isAdded;
    }

    public boolean addEdge(double p, Node startNode, Node endNode) throws NodeNotFoundException {
        if ((!nodes.contains(startNode)) || (!nodes.contains(endNode))) {
            throw new NodeNotFoundException("Node not found!");
        }

        if (p >= Math.random()) {
            Edge e = new Edge(startNode, endNode);
            if (startNode.getEdgeTo(endNode) != null) {
                return false;
            } else {
                startNode.addEdge(e);
                endNode.addEdge(e);
                this.edges.add(e);
                return true;
            }
        }
        return false;
    }

    public void removeEdge (Node startNode, Node endNode) throws NodeNotFoundException {
        if ((!nodes.contains(startNode)) || (!nodes.contains(endNode)))
            throw new NodeNotFoundException("Node not found!");
        edges.remove(this.getEdgeBetween(startNode, endNode));
    }

    public Edge getEdgeBetween(Node startNode, Node endNode) {
        for (Edge edge : edges) {
            if (edge.getStartNode().equals(startNode) && edge.getEndNode().equals(endNode)) {
                return edge;
            }
        }
        return null;
    }

    public Node getNodeByIndex(int index) {
        for (Node n : this.nodes) {
            if (n.getIndex() == index) {
               return n;
            }
        }
        return null;
    }

    public List<Node> getNodes() {
        return nodes;
    }

    public double getNodeAttendanceRate() {
        double sum = 0.0;
        for (Node node : nodes) {
            sum += node.getP();
        }
        return sum / this.nodes.size();
    }

    public List<Edge> getEdges() {
        return edges;
    }

}