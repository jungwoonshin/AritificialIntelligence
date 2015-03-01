/**
 *  
 * CS 440/640 Artificial Intelligence
 * Programming Assignment #2
 * 
 * Masaya Ando, Marika Lee, and Jungwoon Shin
 * February 25, 2014
 * 
 * File: NeuralNet.java
 * 
 * Original Skeleton Code by: Zhiaing Ren (Feb 4, 2012)
 * 
 */

package cs440;

import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;
import java.util.Random;

public class NeuralNet {
	/*
	 * layers: array of the number of nodes in each layer (input and output are
	 * also layers)
	 * 
	 * all indices start from 0
	 */
	public NeuralNet(int[] layers) throws RuntimeException {
		if (layers.length < 2) {
			throw new RuntimeException("The NeuralNet must have at least two layers.");
		}
		m_layers = new ArrayList<List<Node>>(layers.length);

		for (int i = 0; i < layers.length; ++i) {
			List<Node> layer = new ArrayList<Node>(layers[i]);
			for (int k = 0; k < layers[i]; ++k) {
				layer.add(new Node(i, k, false));
			}
			m_layers.add(layer);
		}
	}

	/*
	 * fully connect all the nodes on each layer
	 */
	public void connectAll() {
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
			for (Node node : cur_layer) {
				addThreshold(node, generator.nextDouble());
			}
			pre_layer = cur_layer;
		}
	}

	public void connectTest() {
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

	
	
	// Back Propagation Algorithm Pseudocode
	// Patrick Henry Winston "Artificial Intelligence (3rd Edition)" page 457 
	// 
	// Pick a rate parameter, r 
	// Until Performance is satisfactory
	// 		For each input
	//			Compute the resulting output
	//			Compute Beta for nodes in the output layer using Betaz = dz - oz;
	// 			Compute Beta for nodes in the hidden layer using Betaj = ∑ w j->k * ok * (1 - oj) * Betak
	//			Compute weight changes for all weights using ∆w i->j = r * oi * oj * (1 - oj) * Betaj
	//		Sum up the weight changes for inputs and update weights
	
	public void train(double[][] inputvs, double[][] outputvs, double rate) throws RuntimeException {
		
		// Iterate through input nodes
		for (int i = 0; i < inputvs.length; i++) {
			
			// Compute the resulting output 
			evaluate(inputvs[i]);
			double target = outputvs[i][outputvs[i].length - 1];
			
			// List of nodes in the output layer
			List<Node> outputLayer = m_layers.get(m_layers.size() - 1);
			
			// 1) Compute Beta for output nodes
			for (Node node : outputLayer) {
				// Calculate delta of target node using rate parameter r
				node.calcBetaOutput(target, rate);
			}
			
			// 2) Compute Beta for "hidden layer" nodes
			List<Node> hiddenLayer;
			if (m_layers.size() > 2) {
				// Number of layers excluding input and output layer
				int numLayers = m_layers.size() - 2;
				
				// Iterate through the layers in bottom-up fashion (start from the end)
				for (int j = numLayers; j >= 0; j--) {
					// List of nodes in the one of the hidden layers
					hiddenLayer = m_layers.get(j);
					
					// Iterate through the nodes of the hidden layer
					for (Node node : hiddenLayer) {
						// Calculate deltaw of nodes of hidden layer
						node.calcBetaHidden(rate);
					}
				}
			}
			
			// 3) Compute weight changes for all weights in current layer
			List<Node> currentLayer;
			for (int j = m_layers.size() - 2; j >= 0; j--) {
				currentLayer = m_layers.get(j);
				
				for (Node node : currentLayer) {
					for (Connection c : node.getInputConnectionList()) {
						c.calcDeltaw(rate);
					}
				}
			}
		}
		
		// 4) Add up the weight changes for all inputs and change the weights
		for (int j = m_layers.size() - 1; j >= 0; j--) {
			List<Node> updateLayer = m_layers.get(j);
			
			for (Node node : updateLayer) {				
				for (Connection c : node.getInputConnectionList()) {
					c.updateWeight();
				}
			}
		}
	}	

	
	
	// This method shall change the input and output of each node.
	public double[] evaluate(double[] inputv) throws RuntimeException {
		if (inputv.length != m_layers.get(0).size()) {
			throw new RuntimeException("incompabile inputv");
		}

		Iterator<List<Node>> iter = m_layers.iterator();
		List<Node> layer = iter.next();

		// input layer
		int i = 0;
		for (Node node : layer) {
			node.setOutput(inputv[i]);
			++i;
		}

		while (iter.hasNext()) {
			layer = iter.next();
			calcOutput(layer);
		}

		// copy result
		double[] output = new double[layer.size()];
		i = 0;
		for (Node node : layer) {
			output[i] = node.getOutput();
			++i;
		}

		return output;
	}

	public double error(double[][] inputvs, double[][] outputvs) throws RuntimeException {
		if (inputvs.length != outputvs.length) {
			throw new RuntimeException("inputvs and outputvs are not of the same length");
		}

		double error = 0;

		for (int i = 0; i < inputvs.length; ++i) {
			if (outputvs[i].length != m_layers.get(m_layers.size() - 1).size()) {
				throw new RuntimeException("incompatible outputs");
			}
			double[] results = evaluate(inputvs[i]);
			for (int j = 0; j < results.length; ++j) {
				error += (results[j] - outputvs[i][j]) * (results[j] - outputvs[i][j]);
			}
		}

		error /= inputvs.length;
		error = Math.pow(error, 0.5);

		return error;
	}

	public double errorrate(double[][] inputvs3, double[][] outputvs3, int id) {
		double accu = 0;
		for (int i = 0; i < inputvs3.length; ++i) {
			double[] inputs3 = inputvs3[i];
			double[] results = evaluate(inputs3);
			double target = outputvs3[i][outputvs3[i].length - 1];
			double ret = results[results.length - 1];
			System.out.println("target is " + target + ", ret is " + ret);

			switch (id) {
				default: 
				case 0: 
					if (1.1 - target > 0.5) { // false
						if (ret > 0.5) { // decide to be true
							++accu;
						}
					} else { // true
						if (ret < 0.5) { // decide to be false
							++accu;
						}
					}
					
					break; 
					
				case 1:
					if (1 - target > 0.6) {
						if (ret > 0.333) {
							++accu;
						}
					} else if (1.6 - target > 1) {
						if (ret < 0.333) {
							++accu;
						}
					} else {
						if (ret < 0.9) {
							++accu;
						}
					}
					break;
					
				case 2:
					if (target<0.375) {
						if (ret > 0.375) {
							++accu;
						}
					} else if (target<0.675) {
						if (ret> 0.675 || ret < 0.375) {
							++accu;
						}
					} else if (target<0.850) {
						if (ret > 0.850 || ret < 0.675) {
							++accu;
						}
					} else {
						if (ret < 0.850) {
							++accu;
						}
					}
					break;
					
			}
		}

		double rate = accu / inputvs3.length;
		System.out.println("error rate is " + accu + "/" + inputvs3.length + " = " + rate + "\n");
		return rate;
	}

	private void calcOutput(List<Node> layer) {
		for (Node node : layer) {
			node.calcOutput();
		}
	}

	private List<List<Node>> m_layers;
}
