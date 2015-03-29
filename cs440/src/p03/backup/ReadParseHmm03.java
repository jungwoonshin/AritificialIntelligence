package p03.backup;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.OptionalDouble;


public class ReadParseHmm03 {
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
		int i=0;
		int N=0; //number of states
		int M=0; //number of observation symbols.
		int T=0; //number of length of sequence
		String[] parsedSpaceString=null;
		String[] list_of_vocabs=null;
		double[][] a_matrix =null;
		double[][] b_matrix =null;
		double[] pi_matrix =null;
		@SuppressWarnings("unused")
		String[] parseSpaceAvalues=null;
		@SuppressWarnings("unused")
		double[][] alpha;
		try {
			String sCurrentLine;
			br = new BufferedReader(new FileReader("/Users/jungwoonshin/git/cs440/cs440/src/p03/sentence.hmm"));

			while ((sCurrentLine = br.readLine()) != null) {
				System.out.println("i: " + i);
				System.out.println("sCurrentLine: " + sCurrentLine);
				parsedSpaceString = sCurrentLine.split(" ");
				switch(i){
				case 0:
					N =Integer.parseInt(parsedSpaceString[0]);
					M =Integer.parseInt(parsedSpaceString[1]);
					T = Integer.parseInt(parsedSpaceString[2]);
					a_matrix = new double[N][N];
					b_matrix = new double [N][M];
					pi_matrix = new double[N];
					break;
				case 1:
					break;
				case 2:
					list_of_vocabs = sCurrentLine.split(" ");
					break;
				case 4:
					a_matrix = extract_amatrix_Values(0, a_matrix,sCurrentLine.split(" "));
					break;
				case 5:
					a_matrix = extract_amatrix_Values(1, a_matrix,sCurrentLine.split(" "));
					break;
				case 6:
					a_matrix = extract_amatrix_Values(2, a_matrix,sCurrentLine.split(" "));
					break;
				case 7:
					a_matrix = extract_amatrix_Values(3, a_matrix,sCurrentLine.split(" "));
					break;
				case 9:
					b_matrix = extract_bmatrix_Values(0, b_matrix,parsedSpaceString);
					break;
				case 10:
					b_matrix = extract_bmatrix_Values(1, b_matrix,parsedSpaceString);
					break;
				case 11:
					b_matrix = extract_bmatrix_Values(2, b_matrix,parsedSpaceString);
					break;
				case 12:
					b_matrix = extract_bmatrix_Values(3, b_matrix,parsedSpaceString);
					break;
				case 14:
					pi_matrix[0] = Double.valueOf(parsedSpaceString[0]);
					pi_matrix[1] = Double.valueOf(parsedSpaceString[1]);
					pi_matrix[2] = Double.valueOf(parsedSpaceString[2]);
					pi_matrix[3] = Double.valueOf(parsedSpaceString[3]);
//					pi_matrix[0] = .25;
//					pi_matrix[1] = .25;
//					pi_matrix[2] =  .25;
//					pi_matrix[3] = .25;
					break;
				}

				i++;
			}
			System.out.println("N: " + N);
			System.out.println("M: " + M);
			System.out.println("T: " + T);
			System.out.println("a_matrix: " + Arrays.deepToString(a_matrix));
			System.out.println("b_matrix: " + Arrays.deepToString(b_matrix));
			System.out.println("pi matrix: " + Arrays.toString(pi_matrix));
			System.out.println("==========================================End of Read File==========================================");

			String obsFilename1 = "/Users/jungwoonshin/git/cs440/cs440/src/p03/example1.obs";
			runObsFile(obsFilename1, N, list_of_vocabs, a_matrix, b_matrix, pi_matrix);
			String obsFilename2 = "/Users/jungwoonshin/git/cs440/cs440/src/p03/example2.obs";
			runObsFile(obsFilename2, N, list_of_vocabs, a_matrix, b_matrix, pi_matrix);






		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (br != null)br.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
			br.close();
		}

	}

	private static void runObsFile(String obsFilename, int N, String[] list_of_vocabs,
			double[][] a_matrix, double[][] b_matrix, double[] pi_matrix)
					throws FileNotFoundException, IOException {
		BufferedReader br=new BufferedReader(new FileReader(obsFilename));
		String sCurrentLine =br.readLine() ;
		int dataset_length = Integer.parseInt(sCurrentLine);

		double[][] alpha;
		double[][] sigma;
		double[][] chai;
		int numWords =0;

		for(int dataSetIndex=0; dataSetIndex<dataset_length; dataSetIndex++){
			System.out.println("=============== Start of Forward Algorithm==============================");

			if((sCurrentLine = br.readLine()) != null) {
				numWords = Integer.parseInt(sCurrentLine);
			}
			String words[] =  br.readLine().split(" ");
			System.out.println("dataset_length: " + dataset_length);
			System.out.println("words: " +Arrays.toString(words));
			System.out.println("numWords: " + numWords);

			int obsIndex[] = new int[numWords];
			for(int j=0; j<numWords; j++){ 
				obsIndex[j] =getIndex(words[j], list_of_vocabs);
				System.out.println("obsIndex["+j+"]: " + obsIndex[j] );
			}

			alpha = getAlpha(N, a_matrix, b_matrix, pi_matrix, numWords,
					obsIndex);
			double Answer = 0.0;
			for(int state=0; state<N; state++){ Answer += alpha[numWords-1][state]; }
			System.out.println("Answer: " + Answer);
			
			System.out.println("=============== End of Forward Algorithm==============================");

			
			System.out.println("===============Start of Viterabi Algorithm==============================");

			/* Viterabi Algorithm */
			sigma = new double[numWords][N];
			chai = new double[numWords][N];
			double[] a_ij_times_simga=new double[N];

			// Step 1 & 2: Recursion
			sigma = getSigma(N, a_matrix, b_matrix, pi_matrix,
					sigma, numWords, obsIndex, a_ij_times_simga);
			chai = getChai(N, a_matrix, sigma, chai, numWords, a_ij_times_simga);
			
//			System.out.println("chai: " + Arrays.deepToString(chai));
//			System.out.println("numWords: " + numWords);
//			System.out.println("sigma.length: " + sigma.length);
//			System.out.println("sigma[0].length: " + sigma[0].length);

			// Step 3: termination step
			int q_star[]=new int[numWords];
			q_star = getQStar(sigma, chai, numWords, q_star);
			System.out.println("q_star[i]: " + Arrays.toString(q_star) );
			System.out.println("===============End of Viterabi Algorithm==============================");


		}









		br.close();
	}

	private static int[] getQStar(double[][] sigma, double[][] chai,
			int numWords, int[] q_star) {
		double p_star = Arrays.stream(sigma[numWords-1]).max().getAsDouble();
		for(double e: sigma[numWords-1]){
			if(e==p_star) break;
			q_star[numWords-1]++;
		}
		System.out.println("p_star: " + p_star);
		
		//Path.
		for(int i=numWords-2;i<=0;i--){
			q_star[i] =  (int) chai[i+1][(int)q_star[i+1]];
		}
		return q_star;
	}

	private static double[][] getChai(int N, double[][] a_matrix, double[][] sigma,
			double[][] chai, int numWords, double[] a_ij_times_simga) {
		//t=2 
		for(int t=1;t<numWords;t++){
			for(int j=0;j<N;j++){
				for(int i=1;i<N;i++){
					a_ij_times_simga[i] = sigma[t-1][i]*a_matrix[i][j];
				}
				double highest_value = Arrays.stream(a_ij_times_simga).max().getAsDouble();
				int highest_index=0;
				for(double s:a_ij_times_simga){
					if(s==highest_value) break;
					highest_index++;
				}
				chai[t][j] = highest_index;
			}
		}
		return chai;
	}

	private static double[][] getSigma(int N, double[][] a_matrix,
			double[][] b_matrix, double[] pi_matrix, double[][] sigma,
			int numWords, int[] obsIndex, double[] a_ij_times_simga) {
		//initialization step, t=1
		for(int state=0; state<N;state++){
			sigma[0][state] = pi_matrix[state] * b_matrix[state][obsIndex[0]];
		}
		System.out.println("sigma After initialization: " + Arrays.deepToString(sigma));

		//recursion step

		//t=2
		for(int t=1;t<numWords;t++){
			for(int j=0;j<N;j++){
				a_ij_times_simga= new double[N];
				for(int state=0;state<N;state++){
//						System.out.println("sigma[t-1][state]*a_matrix[state][j]:  " + sigma[t-1][state]*a_matrix[state][j]);
					a_ij_times_simga[state] = sigma[t-1][state]*a_matrix[state][j];
				}
				OptionalDouble highest = Arrays.stream(a_ij_times_simga).max();
				sigma[t][j] = highest.getAsDouble() * b_matrix[j][obsIndex[t]];
			}
		}
		System.out.println("sigma After recursion: " + Arrays.deepToString(sigma));

		return sigma;
	}

	private static double[][] getAlpha(int N, double[][] a_matrix,
			double[][] b_matrix, double[] pi_matrix, int numWords,
			int[] obsIndex) {
		double[][] alpha= new double[numWords][N];
		//Step : t=1
		for(int state=0;state<N;state++){
			//intialize pi prob.
			//initalize step t=1
			alpha[0][state] = pi_matrix[state]*b_matrix[state][obsIndex[0]];
		}
		//Step : t=2,
		for(int t=1; t<numWords;t++){
			for(int l=0; l<N;l++){
				for(int k =0;k<N;k++){
					alpha[t][l] += alpha[t-1][k] * a_matrix[k][l] * b_matrix[l][obsIndex[t]];
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

	private static double[][] extract_amatrix_Values(int i, double[][] a_matrix,
			String[] parseSpaceAvalues) {
		//		System.out.println("Double.valueOf(parseSpaceAvalues[0]):  " +Double.valueOf(parseSpaceAvalues[0]));
		//		System.out.println("Double.valueOf(parseSpaceAvalues[1]):  " +Double.valueOf(parseSpaceAvalues[1]));
		//		System.out.println("Double.valueOf(parseSpaceAvalues[2]):  " +Double.valueOf(parseSpaceAvalues[2]));
		//		System.out.println("Double.valueOf(parseSpaceAvalues[3]):  " +Double.valueOf(parseSpaceAvalues[3]));
		//
		//		a_matrix[i][0] = Double.valueOf(parseSpaceAvalues[0]);
		//		a_matrix[i][1] = Double.valueOf(parseSpaceAvalues[1]);
		//		a_matrix[i][2] = Double.valueOf(parseSpaceAvalues[2]);
		//		a_matrix[i][3] = Double.valueOf(parseSpaceAvalues[3]);
		for(int k=0; k<a_matrix[0].length;k++){
			a_matrix[i][k] = Double.valueOf(parseSpaceAvalues[k]);
		}
		return a_matrix;
	}

	private static double[][] extract_bmatrix_Values(int i, double[][] b_matrix,
			String[] parseSpaceAvalues) {
		for(int k=0;k<b_matrix[0].length;k++){
			//			System.out.println("Double.valueOf(parseSpaceAvalues[k]):  " +Double.valueOf(parseSpaceAvalues[k]));
			b_matrix[i][k] = Double.valueOf(parseSpaceAvalues[k]);
		}
		return b_matrix;
	}

}
