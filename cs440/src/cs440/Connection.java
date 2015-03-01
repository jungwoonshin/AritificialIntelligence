/**
 * @author Zhiqiang Ren 
 * date: Feb. 4th. 2012
 * 
 */

package cs440;


public class Connection {
    
    public Connection(Node from, Node to, double weight) {
        m_from = from;
        m_to = to;
        m_weight = weight;
    }
    
    public Node getFromNode() {
        return m_from;
    }
    
    public Node getToNode() {
        return m_to;
    }
    
    public double getWeight() {
        return m_weight;
    }
    
    public void updateWeight() {
    	m_weight += m_deltaw;
    	m_deltaw = 0;
    }
    
    public void calcDeltaw(double rate) {
    	m_deltaw += rate * m_from.getOutput() * m_to.getOutput() * (1 - m_to.getOutput()) * m_to.getBeta();
    }
        
    private double m_weight;
    private double m_deltaw;

    private Node m_from;
    private Node m_to;

}
