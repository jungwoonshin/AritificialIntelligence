/**
 * @author Zhiqiang Ren 
 * date: Feb. 4th. 2012
 * 
 */
package aipackage;

import java.io.FileNotFoundException;
import java.util.List;
import java.util.Random;

/**
 * @author Zhiqiang Ren
 * 
 */
public class NeuralNetLearner {
	/**
	 * @param args
	 * @throws FileNotFoundException
	 */
	public static void main(String[] args) throws FileNotFoundException {
		int[] layers = { 6, 2, 1 }; // three layers
		NeuralNet net = new NeuralNet(layers);


		net.connectTest();

		//First layer
		List<Node> node1 = net.m_layers.get(0);

		//Witin First layer, get all the nodes.
		for(int i=0; i<node1.size(); i++){

			//accessing nodes
			Node layer0_nodeI = node1.get(i);
			
			System.out.println(layer0_nodeI.getInputConnection().toString());
			System.out.println(layer0_nodeI.getOutputConnection(0).toString());

			//System.out.println(output);
		}




		double[][] inputvs = { { 1, 1, 0, 0, 0, 0 }, { 1, 0, 1, 0, 0, 0 },
				{ 1, 0, 0, 1, 0, 0 }, { 1, 0, 0, 0, 1, 0 },
				{ 1, 0, 0, 0, 0, 1 }, { 0, 1, 1, 0, 0, 0 },
				{ 0, 1, 0, 1, 0, 0 }, { 0, 1, 0, 0, 1, 0 },
				{ 0, 1, 0, 0, 0, 1 }, { 0, 0, 1, 1, 0, 0 },
				{ 0, 0, 1, 0, 1, 0 }, { 0, 0, 1, 0, 0, 1 },
				{ 0, 0, 0, 1, 1, 0 }, { 0, 0, 0, 1, 0, 1 },
				{ 0, 0, 0, 0, 1, 1 } };

		double[][] outputvs = { { 0 }, { 0 }, { 1 }, { 1 }, { 1 }, { 0 },
				{ 1 }, { 1 }, { 1 }, { 1 }, { 1 }, { 1 }, { 0 }, { 0 }, { 0 } };

		//        for (int n = 0; n < 300; ++n) {
		//            net.train(inputvs, outputvs, 1);
		//        }
		//
		//        net.errorrate(inputvs, outputvs);
		//        System.out.println("============================");
		//
		//        int[] layers2 = { 2, 2 }; // two layers
		//        NeuralNet net2 = new NeuralNet(layers2);
		//        net2.connectAll();
		//
		//        double[][] inputvs2 = { { 0, 0 }, { 0, 1 }, { 1, 1 }, { 1, 0 } };
		//        double[][] outputvs2 = { { 0, 0 }, { 0, 1 }, { 1, 1 }, { 0, 1 } };
		//
		//        for (int n = 0; n < 100; ++n) {
		//            net2.train(inputvs2, outputvs2, 1);
		//        }
		//        
		//        net2.errorrate(inputvs2, outputvs2);
		//        System.out.println("============================");
		//
		//        DataProcessor data = new DataProcessor("crx.data.training", 0 );
		//        int[] layers3 = { 15, 30, 1 }; // two layers
		//        NeuralNet net3 = new NeuralNet(layers3);
		//        net3.connectAll();
		//        
		//        double[][] inputvs3 = data.m_inputvs;
		//        double[][] outputvs3 = data.m_outputvs;
		//
		//        for (int n = 0; n < 10; ++n) {
		//            net3.train(inputvs3, outputvs3, 3);
		//
		//            double error = net3.error(inputvs3, outputvs3);
		//            System.out.println("error is " + error);
		//        }
		//        
		//        net3.errorrate(inputvs3, outputvs3);

		return;
	}

}
