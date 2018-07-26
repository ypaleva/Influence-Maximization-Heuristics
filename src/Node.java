import java.text.DecimalFormat;
import java.util.*;

public class Node {

    private int index;
    private double p;
    private double degree;
    private double averageCoverage;
    private double rewiringProb;
    private boolean invited;
    private boolean attended;
    private boolean covered;
    private List<Edge> edgesFrom;
    private List<Edge> edgesTo;
    private int uniqueNodesCovered;

    public Node(int index, double p, double degree, double rewiringProb, boolean invited, boolean attended, int uniqueNodesCovered) {
        this.index = index;
        this.p = p;
        this.degree = degree;
        this.rewiringProb = rewiringProb;
        this.invited = invited;
        this.attended = attended;
        this.edgesFrom = new ArrayList<>();
        this.edgesTo = new ArrayList<>();
        this.uniqueNodesCovered = uniqueNodesCovered;
    }

    public void addEdgeTo(Node nodeTo) {
        Edge edgeTo = new Edge(this, nodeTo);
        this.edgesTo.add(edgeTo);
    }

    public void addEdgeFrom(Node nodeFrom) {
        Edge edgeFrom = new Edge(this, nodeFrom);
        this.edgesFrom.add(edgeFrom);
    }

    public boolean addEdge(Edge edge) {
        if (edge.getStartNode().equals(this)) {
            this.edgesFrom.add(edge);
        } else if (edge.getEndNode().equals(this)) {
            this.edgesTo.add(edge);
        } else {
            return false;
        }
        return false;
    }

    public int getNumEdgesTo() {
        return this.edgesTo.size();
    }

    public Edge getEdgeTo(Node nodeTo) {
        for (Edge e : this.edgesFrom) {
            if (e.getEndNode().equals(nodeTo)) {
                return e;
            }
        }
        return null;
    }

    public int getIndex() {
        return this.index;
    }

    public double getP() {
        return p;
    }

    public boolean isInvited() {
        return invited;
    }

    public double getDegree() {
        return degree;
    }

    public void setDegree(double degree) {
        this.degree = degree;
    }

    public void setIndex(int index) {
        this.index = index;
    }

    public void setP(double p) {
        this.p = p;
    }

    public void setInvited(boolean invited) {
        this.invited = invited;
    }

    public boolean isAttended() {
        return attended;
    }

    public void setAttended(boolean attended) {
        this.attended = attended;
    }

    public List<Edge> getEdgesFrom() {
        return edgesFrom;
    }

    public void setEdgesFrom(List<Edge> edgesFrom) {
        this.edgesFrom = edgesFrom;
    }

    public void setEdgesTo(List<Edge> edgesTo) {
        this.edgesTo = edgesTo;
    }

    public List<Node> getNodesReachable() {
        List<Node> nodesConnected = new ArrayList<>();
        for (Edge edge : this.edgesFrom) {
            nodesConnected.add(edge.getEndNode());
        }
        return nodesConnected;
    }

    public List<Integer> getNodesConnected() {
        List<Integer> nodesConnected = new ArrayList<>();
        for (Edge edge : this.edgesFrom) {
            nodesConnected.add(edge.getEndNode().getIndex());
        }
        return nodesConnected;
    }

    public double getRewiringProb() {
        return rewiringProb;
    }

    public double getBestAverageCoverage() {
        return this.getAverageCoverage() * this.getP() ;
    }


    public double getAverageCoverage() {
        return averageCoverage ;
    }

    public void setAverageCoverage(double averageCoverage) {
        this.averageCoverage = averageCoverage;
    }

    public void setRewiringProb(double rewiringProb) {
        this.rewiringProb = rewiringProb;
    }

    public int getNeighboursSize() {
        return getNodesReachable().size();
    }

    public List<Edge> getEdgesTo() {
        return edgesTo;
    }

    public Integer getUniqueNodesCovered() {
        return this.uniqueNodesCovered;
    }

    public void setUniqueNodesCovered(int uniqueNodesCovered) {
        this.uniqueNodesCovered = uniqueNodesCovered;
    }

    public double getAddedValue() {
        return this.getUniqueNodesCovered() * this.getP();
    }

    public boolean isCovered() {
        return covered;
    }

    public void setCovered(boolean covered) {
        this.covered = covered;
    }

    public double distribution() {
        return this.getP() * this.getNeighboursSize();
    }

    private static double roundThreeDecimals(double d) {
        DecimalFormat twoDForm = new DecimalFormat("#.###");
        return Double.valueOf(twoDForm.format(d));
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("Node{");
        sb.append("i=").append(index);
        sb.append(", p=").append(roundThreeDecimals(p));
        sb.append(", n=").append(this.getNodesReachable().size());
        sb.append('}');
        return sb.toString();
    }


    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Node node = (Node) o;

        if (getIndex() != node.getIndex()) return false;
        if (Double.compare(node.getP(), getP()) != 0) return false;
        if (isInvited() != node.isInvited()) return false;
        if (isAttended() != node.isAttended()) return false;
        if (getUniqueNodesCovered() != node.getUniqueNodesCovered()) return false;
        if (getEdgesFrom() != null ? !getEdgesFrom().equals(node.getEdgesFrom()) : node.getEdgesFrom() != null)
            return false;
        return getEdgesTo() != null ? getEdgesTo().equals(node.getEdgesTo()) : node.getEdgesTo() == null;
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        result = getIndex();
        temp = Double.doubleToLongBits(getP());
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        result = 31 * result + (isInvited() ? 1 : 0);
        result = 31 * result + (isAttended() ? 1 : 0);
        result = 31 * result + (getEdgesFrom() != null ? getEdgesFrom().hashCode() : 0);
        result = 31 * result + (getEdgesTo() != null ? getEdgesTo().hashCode() : 0);
        result = 31 * result + getUniqueNodesCovered();
        return result;
    }
}