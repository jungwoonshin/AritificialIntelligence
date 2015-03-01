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

package p2;

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
	
		System.out.println("============= Test Neural Net 1 ===============");

		int numTraining1 = 3; 

		int[] layers = { 6, 2, 1 }; // three layers
		NeuralNet net = new NeuralNet(layers);
		
		
		allowAccesstoNeuralNetFromNode(layers, net);
		
		
		
//		net.connectTest();
		net.connectAll();
		double[][] inputvs = { { 1, 1, 0, 0, 0, 0 }, { 1, 0, 1, 0, 0, 0 },
				{ 1, 0, 0, 1, 0, 0 }, { 1, 0, 0, 0, 1, 0 },
				{ 1, 0, 0, 0, 0, 1 }, { 0, 1, 1, 0, 0, 0 },
				{ 0, 1, 0, 1, 0, 0 }, { 0, 1, 0, 0, 1, 0 },
				{ 0, 1, 0, 0, 0, 1 }, { 0, 0, 1, 1, 0, 0 },
				{ 0, 0, 1, 0, 1, 0 }, { 0, 0, 1, 0, 0, 1 },
				{ 0, 0, 0, 1, 1, 0 }, { 0, 0, 0, 1, 0, 1 },
				{ 0, 0, 0, 0, 1, 1 } };

		double[][] outputvs = { { 0 }, { 0 }, { 1 }, { 1 }, { 1 }, 
				{ 0 }, { 1 }, { 1 }, { 1 }, { 1 }, 
				{ 1 }, { 1 }, { 0 }, { 0 }, { 0 } };

		
		printAllWeightValues(net);
		
		
//		net.train(inputvs, outputvs, 1);
		
		
		for (int n = 0; n < numTraining1; ++n) {
			net.train2(inputvs, outputvs, 1);
		}

		System.out.println(net.errorrate(inputvs, outputvs, 0));
		  //learning test.


		System.out.println("============= Test Neural Net 2 ===============");

		int numTraining2 = 350;

		int[] layers2 = { 2, 2 }; // two layers
		NeuralNet net2 = new NeuralNet(layers2);
		net2.connectAll();

		double[][] inputvs2  = { { 0, 0 }, { 0, 1 }, { 1, 1 }, { 1, 0 } };
		double[][] outputvs2 = { { 0, 0 }, { 0, 1 }, { 1, 1 }, { 0, 1 } };

		for (int n = 0; n < numTraining2; ++n) {
			net2.train(inputvs2, outputvs2, 1);
		}

		net2.errorrate(inputvs2, outputvs2, 0);




		System.out.println("============= Credit Data Training Data===============");

		int numTrainingCredit = 500;

		DataProcessor dataCredit = new DataProcessor("crx.data.training", 0);
		int[] layers3 = { 15, 30, 1 }; // three layers
		NeuralNet net3 = new NeuralNet(layers3);
		net3.connectAll();

		double[][] inputvs3 = dataCredit.m_inputvs;
		double[][] outputvs3 = dataCredit.m_outputvs;


		for (int n = 0; n < numTrainingCredit; ++n) {
			net3.train(inputvs3, outputvs3, .5);
			// double error = net3.error(inputvs3, outputvs3);
			// System.out.println("error is " + error);
		}

		System.out.println(net3.errorrate(inputvs3, outputvs3, 0));





		System.out.println("============= Credit Data Testing Data===============");

		dataCredit = new DataProcessor("crx.data.testing", 0);
		net3.connectAll();
		inputvs3 = dataCredit.m_inputvs;
		outputvs3 = dataCredit.m_outputvs;

		System.out.println(net3.errorrate(inputvs3, outputvs3, 0));




		System.out.println("============= Lens Data Training Data===============");

		int numTraining4 = 2000;

		DataProcessor dataLens = new DataProcessor("lenses.training", 1);
		int[] layers4 = { 4, 30, 1 };
		NeuralNet net4 = new NeuralNet(layers4);
		net4.connectAll();

		double[][] inputvs4 = dataLens.m_inputvs;
		double[][] outputvs4 = dataLens.m_outputvs;

		for (int i = 0; i < numTraining4; i++) {
			net4.train(inputvs4, outputvs4, 1);
		}

		System.out.println(net4.errorrate(inputvs4, outputvs4, 1));

		System.out.println("============= Lens Data Training Data===============");

		dataLens = new DataProcessor("lenses.testing", 1);

		net4.connectAll();

		inputvs4 = dataLens.m_inputvs;
		outputvs4 = dataLens.m_outputvs;


		System.out.println(net4.errorrate(inputvs4, outputvs4, 1));
	
		
		
		/* Bubil Data Set.
		for(int k=2;k<31;k++){
			for(int j=0;j<20;j++){

				//		int k=4; int j=9;

//				System.out.println("============= BUBIL Data Training Data===============");


				int numTraining5 = 2000;

				DataProcessor dataBubils = new DataProcessor("BUBIL.training", 2);
				int[] layers5 = { 4, 30,1 };
				layers5[1] = k;
				NeuralNet net5 = new NeuralNet(layers5);
				net5.connectAll();

				double[][] inputvs5 = dataBubils.m_inputvs;
				double[][] outputvs5 = dataBubils.m_outputvs;

				//System.out.println("===k: " + k + "===j : " + j + "==============");

				for (int i = 0; i < numTraining5; i++) {
					net5.train(inputvs5, outputvs5, j);
				}
				double rate1 = net5.errorrate(inputvs5, outputvs5, 2);

				//				System.out.println("============= BUBIL Data testing Data===============");

				dataBubils = new DataProcessor("BUBIL.testing", 2);

				net5.connectAll();

				inputvs5 = dataBubils.m_inputvs;
				outputvs5 = dataBubils.m_outputvs;

				double rate2 = net5.errorrate(inputvs5, outputvs5, 2);
				if (rate2 <0.5)
					System.out.println("rate1: " +rate1 +", rate2: " + rate2 +", k:" + k +", j:" + j);;
					//		System.out.println("rate1: " + rate1 + ", rate2: " + rate2);

			}//inner loop
		} 
		 */


		return;
	}

	private static void allowAccesstoNeuralNetFromNode(int[] layers,
			NeuralNet net) {
		for (int i = 0; i < layers.length; ++i) {
			List<Node> layer = net.m_layers.get(i);
			for (int k = 0; k < layers[i]; ++k) {
				layer.get(k).setNeuralNet(net);
			}
		}
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
					System.out.println("layer:" + i + " , Node:" + j+ " to Node:" + k + "'s weight value: " + outCN.getWeight());
					k++;
				}
				j++;
				System.out.println();
			}
		}
		
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
					System.out.println("layer:" + i + " , Node:" + k+ " to Node:" + j + "'s weight value: " + inCN.getWeight());
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
					System.out.println("layer:" + i + " , Node:" + k+ " to Node:" + j + "'s weight value: " + inCN.getWeight());
					
				}
				System.out.println();
				
				
			}
			System.out.println();
		
		}
		
		
		
		
		
	}

}
