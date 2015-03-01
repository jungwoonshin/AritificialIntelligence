/**
 * @author Zhiqiang Ren 
 * date: Feb. 4th. 2012
 * 
 */

package aipackage;


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
    
    private double m_weight;
    private double m_deltaw;

    private Node m_from;
    private Node m_to;
	@Override
	public String toString() {
		return "Connection [m_from=" + m_from + ", m_to=" + m_to
				+ ", m_weight=" + m_weight + ", m_deltaw=" + m_deltaw + "]";
	}
	
    
}
