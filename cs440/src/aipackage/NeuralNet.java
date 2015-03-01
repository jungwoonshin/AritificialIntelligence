/**
 * @author Zhiqiang Ren 
 * date: Feb. 4th. 2012
 * 
 */
package aipackage;

import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;
import java.util.Random;

public class NeuralNet {
    /*
     * layers: array of the number of nodes in each layer (input and output are also layers)
     * 
     * all indices start from 0
     */

    public List<List<Node> > m_layers;

    
    public NeuralNet(int [] layers) throws RuntimeException
    {
        if (layers.length < 2)
        {
            throw new RuntimeException("The NeuralNet must have at least two layers.");
        }
        m_layers = new ArrayList<List<Node>>(layers.length);
        
        for (int i = 0; i < layers.length; ++i)
        {
            List<Node> layer = new ArrayList<Node>(layers[i]);
            for (int k = 0; k < layers[i]; ++k)
            {
                layer.add(new Node(i, k, false));
            }
            m_layers.add(layer);
        }
    }
    
    /*
     * fully connect all the nodes on each layer
     */
    public void connectAll()
    {
        Random generator = new Random();
        Iterator<List<Node>> iter = m_layers.iterator();
        
        List<Node> pre_layer = iter.next();
        while (iter.hasNext()) {
            List<Node> cur_layer = iter.next();

            for (int i = 0; i < pre_layer.size(); ++i) {
                for (int j = 0; j < cur_layer.size(); ++j) {
                    addConnection(pre_layer, i, cur_layer, j, generator.nextDouble());
                }
            }
            for (Node node: cur_layer) {
                addThreshold(node, generator.nextDouble());
            }
            pre_layer = cur_layer;
        }
    }
    
    public void connectTest()
    {
        addConnection(0, 0, 1, 0, 0.2);
        addConnection(0, 1, 1, 0, 0.3);
        addConnection(0, 2, 1, 0, 0.4);
        addConnection(0, 3, 1, 1, 0.6);
        addConnection(0, 4, 1, 1, 0.7);
        addConnection(0, 5, 1, 1, 0.8);
        
        addConnection(1, 0, 2, 0, 1.0);
        addConnection(1, 1, 2, 0, 1.1);
        
        addThreshold(1, 0, 0.1);
        addThreshold(1, 1, 0.5);
        addThreshold(2, 0, 0.9);
    }
    
    void addConnection(int from_layer, int from_pos, int to_layer, int to_pos, double weight) {
        List<Node> layer_f = m_layers.get(from_layer);
        List<Node> layer_t = m_layers.get(to_layer);
        addConnection(layer_f, from_pos, layer_t, to_pos, weight);
    }
    
    void addConnection(List<Node> layer_f, int from_pos, List<Node> layer_t, int to_pos, double weight) {
        Node from_node = layer_f.get(from_pos);
        Node to_node = layer_t.get(to_pos);
        addConnection(from_node, to_node, weight);
    }
    
    void addConnection(Node from_node, Node to_node, double weight) {
        Connection con = new Connection(from_node, to_node, weight);
        from_node.addOutputConnection(con);
        to_node.addInputConnection(con);
    }
    
    // add a threshold to certain node
    void addThreshold(int layer, int pos, double weight) {
        List<Node> layer_i = m_layers.get(layer);
        Node node = layer_i.get(pos);
        addThreshold(node, weight);
    }
    
    void addThreshold(Node node, double weight) {        
        Node thrd = new Node(node.getLayer(), node.getPos(), true);
        thrd.setOutput(-1);
        
        Connection con = new Connection(thrd, node, weight);
        
        thrd.addOutputConnection(con);
        node.addInputConnection(con);
    }
    
    // r: rate parameter
    public void train(double [][] inputvs, double [][] outputvs, double r) throws RuntimeException
    {
    	
    	
        // todo add your code here
    }
        
    // This method shall change the input and output of each node.
    public double [] evaluate(double [] inputv) throws RuntimeException {
        if (inputv.length != m_layers.get(0).size()) {
            throw new RuntimeException("incompabile inputv");
        }
        
        Iterator<List<Node>> iter = m_layers.iterator();
        List<Node> layer = iter.next();
        
        // input layer
        int i = 0;
        for (Node node: layer) {
            node.setOutput(inputv[i]);
            ++i;
        }
        
        while (iter.hasNext()) {
            layer = iter.next();
            calcOutput(layer);
        }

        // copy result
        double [] output = new double [layer.size()];
        i = 0;
        for (Node node: layer) {
            output[i] = node.getOutput();
            ++i;
        }
        
        return output;
    }
    
    public double error(double [][] inputvs, double [][] outputvs) throws RuntimeException
    {
        if (inputvs.length != outputvs.length)
        {
            throw new RuntimeException("inputvs and outputvs are not of the same length");
        }
        
        double error = 0;
        
        for (int i = 0; i < inputvs.length; ++i) {
            if (outputvs[i].length != m_layers.get(m_layers.size() - 1).size()) {
                throw new RuntimeException("incompatible outputs");
            }
            double [] results = evaluate(inputvs[i]);
            for (int j = 0; j < results.length; ++j) {
                error += (results[j] - outputvs[i][j]) * (results[j] - outputvs[i][j]);
            }
        }
        
        error /= inputvs.length;
        error = Math.pow(error, 0.5);
        
        return error;
    }
        
    public double errorrate(double [][]inputvs3, double [][]outputvs3) {
        double accu = 0;
        for (int i = 0; i < inputvs3.length; ++i) {
            double [] inputs3 = inputvs3[i];
            double [] results = evaluate(inputs3);
            double target = outputvs3[i][outputvs3[i].length - 1];
            double ret = results[results.length - 1];
//            System.out.println("target is " + target + ", ret is " + ret);
            
            if (1.1 - target > 0.5) {  // false
                if (ret > 0.5) {  // decide to be true
                    ++accu;
                }
            } else {  // true
                if (ret < 0.5) {  // decide to be false
                    ++accu;
                }
            }
        }
        
        double rate = accu / inputvs3.length;
        System.out.println("error rate is " + accu + "/" + inputvs3.length + " = " + rate);
        return rate;
      }
    
    private void calcOutput(List<Node> layer) {
        for (Node node: layer) {
            node.calcOutput();
        }           
    }
    
}
