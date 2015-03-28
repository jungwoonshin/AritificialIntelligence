/**
 *  
 * CS 440/640 Artificial Intelligence
 * Programming Assignment #2
 * 
 * Masaya Ando, Marika Lee, and Jungwoon Shin
 * February 25, 2014
 * 
 * File: NeuralNetLearner.java
 * 
 * Original Skeleton Code by: Zhiaing Ren (Feb 4, 2012)
 * 
 */

package new01;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


public class NeuralNetLearner01 {
	/**
	 * @param args
	 * @throws FileNotFoundException
	 */
	public static void main(String[] args) throws FileNotFoundException {
	
		
		testNeuralNet2();

		return;
	}

	static void allowAccesstoNeuralNetFromNode(int[] layers,
			NeuralNet net) {
		for (int i = 0; i < layers.length; ++i) {
			List<Node> layer = net.m_layers.get(i);
			for (int k = 0; k < layers[i]; ++k) {
				layer.get(k).setNeuralNet(net);
			}
		}
	}

	private static void testNeuralNet2() {
		System.out.println("============= Test Neural Net 2 ===============");

		int numTraining2 = 350;

		int[] layers2 = { 2, 2,1 }; // two layers
		NeuralNet net2 = new NeuralNet(layers2);
		net2.connectAll();

		double[][] inputvs2  = { { 0, 0 }, { 0, 1 }, { 1, 1 }, { 1, 0 } };
		double[][] outputvs2 = { { 0, 0 }, { 0, 1 }, { 1, 1 }, { 0, 1 } };

		
		NeuralNetLearner01.printAllWeightValues(net2);

		System.out.println("Beginning to train");
		net2.train2(inputvs2,outputvs2,1);
		NeuralNetLearner01.printAllWeightValues(net2);
		/*
		for (int n = 0; n < numTraining2; ++n) {
			net2.train(inputvs2, outputvs2, 1);
		}
		*/

//		System.out.println(net2.errorrate(inputvs2, outputvs2, 0));
	}
	static double printTwoDigits(double a){
		return ((int)(a*100))/100.0;
		
	}
	
	static void printAllWeightValues(NeuralNet net) {
		System.out.println("Output connection lists");
		for (int i=0; i<net.m_layers.size(); i++){
			List<Node> layer = net.m_layers.get(i);
			int j=0;
			for (Node node: layer) {		
				Iterator<Connection> iter = node.getOutputConnectionList().iterator();
				Connection outCN=null;
				double weight; int k=0;
				while(iter.hasNext()){
					outCN = iter.next();
					weight = outCN.getWeight();
					System.out.println("layer:(" + i + ") , Node:(" + j+ ", " + k + "),  weight: " + printTwoDigits(outCN.getWeight())) ;
					System.out.println("Node output: " + layer.get(j));
					k++;
				}
				j++;
				System.out.println();
			}
		}
		
//		System.out.println("double test: " + printTwoDigits(12.3344353));
		
		System.out.println("Input connection lists1");
		/*
		for (int i=1; i<net.m_layers.size();i++){
			List<Node> layer = net.m_layers.get(i);
			int j=0;
			for (Node node: layer) {		
				Iterator<Connection> iter = node.getInputConnectionList().iterator();
				Connection inCN=null;
				double weight; int k=0;
				while(iter.hasNext()){
					inCN = iter.next();
					weight = inCN.getWeight();
					System.out.println("layer:" + i + " , Node(" + k+ " to Node:" + j + "'s weight value: " + inCN.getWeight());
					k++;
				}
				j++;
				System.out.println();
			}
		}
		*/
		
		
		//Notice that this doesn't use up the threshold node by using size-1.
		for (int i=1; i<net.m_layers.size();i++){
			List<Node> layer = net.m_layers.get(i);
			for(int j=0; j<layer.size();j++){
				Node node=layer.get(j);
				List<Connection> connections = node.getInputConnectionList();
				for(int k=0;k<connections.size()-1;k++){
					Connection inCN = connections.get(k);
					System.out.println("layer:(" + i + ") , Node:(" + k+ ", " + j + "),  weight: " + printTwoDigits(inCN.getWeight()));
					
				}
				System.out.println();
				
				
			}
			System.out.println();
		
		}
		
		
		
		
		
	}

}
