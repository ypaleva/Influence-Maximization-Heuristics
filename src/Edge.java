public class Edge {

    private Node startNode;
    private Node endNode;

    public Edge(Node startNode, Node endNode) {
        this.startNode = startNode;
        this.endNode = endNode;
    }

    public Node getStartNode() {
        return startNode;
    }

    public Node getEndNode() {
        return endNode;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("Edge{");
        sb.append("startNode = ").append(startNode.getIndex());
        sb.append(", endNode = ").append(endNode.getIndex());
        sb.append('}');
        return sb.toString();
    }
}