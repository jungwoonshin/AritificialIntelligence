
package p03.backup;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

public class testWrittenHomework01 {
	public static void main(String[] args) throws Exception {
		BufferedReader br = null;
		/*
		The format of an .hmm file is as follows:  
			The first line contains integers N (number of states), M (number of observation symbols), and T (number of time steps or length of oberservation sequences).  
			The second contains four strings 
			SUBJECT AUXILIARY PREDICATE OBJECT ,
			which refer to four basic English syntactic structures. Each is used to name an individual HMM state.
			The third line contains strings
			kids robots do can play eat chess food, 
			that provide the vocabulary to be used in the observation sentences. 
			Then comes a line with the text "a:", followed by the matrix a.   The matrix b and vector pi are similarly represented.  The matrix and vector elements are floating-point numbers less than or equal to 1.0.
		 */




		try 
		{

			int i=0;
			int N=0; //number of states
			int M=0; //number of observation symbols.
			int T=0; //number of length of sequence
			String[] parsedSpaceString=null;
			String[] states=null;
			String[] list_of_vocabs=null;
			double[][] a_matrix ={{1./2,1./4,1./4},{0,1./2,1./2},{0.,0.,1.}};
			double[][] b_matrix ={{3./4,1./4},{1./4,3./4},{1./2,1./2}};
			double[] pi_matrix ={3./4,0.,1./4};
			String[] parseSpaceAvalues=null;
			double[][] alpha = null;


		
			
			
			N =3;
			M =2;

			System.out.println("N: " + N);
			System.out.println("M: " + M);
			System.out.println("T: " + T);
			System.out.println("a_matrix: " + Arrays.deepToString(a_matrix));
			System.out.println("b_matrix: " + Arrays.deepToString(b_matrix));
			System.out.println("==========================================End of Read File==========================================");

			String obsFilename1 = "/Users/jungwoonshin/git/cs440/cs440/src/p03/example1.obs";
			runObsFile(obsFilename1, N, list_of_vocabs, a_matrix, b_matrix, pi_matrix);
			//			String obsFilename2 = "/Users/jungwoonshin/git/cs440/cs440/src/p03/example2.obs";
			//			runObsFile(obsFilename2, N, list_of_vocabs, a_matrix, b_matrix, pi_matrix);






		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (br != null)br.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}

	}

	private static void runObsFile(String obsFilename, int N, String[] list_of_vocabs,
			double[][] a_matrix, double[][] b_matrix, double[] pi_matrix)
					throws FileNotFoundException, IOException {
		BufferedReader br;
		String sCurrentLine;
		double[][] alpha;


		int t = 2;

		System.out.println("a_matrix: " + Arrays.deepToString(a_matrix));
		System.out.println("b_matrix: " + Arrays.deepToString(b_matrix));

		System.out.println("t: " + t);
		int obsIndex[] = {0,1};


		alpha = forwardProcedure(N, a_matrix, b_matrix, pi_matrix, t,
				obsIndex);
		double seqProb = 0.0;
		for(int state=0; state<N; state++){ seqProb += alpha[t-1][state]; }
		System.out.println("seqProb: " + seqProb);
		System.out.println("==============");

	}
	private static double[][] forwardProcedure(int N, double[][] a_matrix,
			double[][] b_matrix, double[] pi_matrix, int numWords,
			int[] obsIndex) {
		double[][] alpha;
		alpha = new double[numWords+1][N];
		//Step : t=1
		for(int state=0;state<N;state++){
			//				System.out.println("obsIndex[1]-1: "+ (obsIndex[1]-1));

			//intialize pi prob.
			//				System.out.println(" pi_matrix[state]: " +  pi_matrix[state]);
			alpha[0][state] = pi_matrix[state];

			//initalize step t=1
			alpha[1][state] = pi_matrix[state]*b_matrix[state][0];
		}
		//Step : t=2,
		for(int t=2; t<numWords+1;t++){
			for(int l=0; l<N;l++){
				for(int k =0;k<N;k++){
					alpha[t][l] += alpha[t-1][k] * a_matrix[k][l] * b_matrix[l][1];
				}
			}
		}

		System.out.println("alpha: "+ Arrays.deepToString(alpha));
		return alpha;
	}
	public static int getIndex(String word, String vocab[]){
		int i = 0;
		while(i < vocab.length && !vocab[i].equalsIgnoreCase(word) ) {
			i++; 
		}
		return i;
	}

}
