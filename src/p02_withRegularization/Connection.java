/**
 * @author Zhiqiang Ren 
 * date: Feb. 4th. 2012
 * 
 */

package p02_withRegularization;


public class Connection {
	  
    private double m_weight;
    private double m_deltaw;

    private Node m_from;
    private Node m_to;
    
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
    
    // i is the index for previous layer
    // j is the index for next layer
    public void calcDeltaw(double rate) {
//    	System.out.println(" m_to.getBeta(): " +  m_to.getBeta());
    	// R * O_I * O_J * (1-O_J) * Beta value of output layer's node
    	m_deltaw += rate * m_from.getOutput() * m_to.getOutput() 
    			* (1 - m_to.getOutput()) * m_to.getBeta()
    			+1/15.0*m_weight;
//    	m_deltaw += rate * m_from.getOutput() * m_to.getOutput() 
//    			* (1 - m_to.getOutput()) * m_to.getBeta();
    			
    	
    }
    


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
}
