package p03;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.OptionalDouble;
import java.util.Random;

public class testWrittenHomework02 {
	private static final double DENOMINATOR_CONSTANT = 10e-70;

	public static void main(String[] args) throws Exception {
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

	}

	private static void runObsFile(String obsFilename, int N, String[] list_of_vocabs,
			double[][] a_matrix, double[][] b_matrix, double[] pi_matrix)
					throws FileNotFoundException, IOException {
		BufferedReader br;
		String sCurrentLine;
		double[][] alpha;
		double[][] beta;
		int numWords = 2;
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
		
		
		
		alpha = getAlpha(N, a_matrix, b_matrix, pi_matrix, numWords,
				obsIndex);
		beta = getBeta(N, a_matrix, b_matrix, numWords, obsIndex);

		System.out.println("alpha: "+ Arrays.deepToString(alpha));
		System.out.println("beta: "+ Arrays.deepToString(beta));

		double Answer = 0.0;
		for(int state=0; state<N; state++){ Answer += alpha[numWords-1][state]; }
//		System.out.println("sum(alpha(T)): The probablity of given obervation data set [O1...OT]="+Arrays.toString(words)+"=" + Answer);
		System.out.println("=============== End of Forward Algorithm==============================");
		double[][] sigma;
		double[][] chai;
		double[][][] xi;
		double[][] gamma;
		int M =2;
		System.out.println("\n\n===============Start of Viterabi Algorithm==============================");
		// Viterabi Algorithm 
		sigma = new double[numWords][N];
		chai = new double[numWords][N];
		double[] a_ij_times_simga=new double[N];

		// Step 1 & 2: Recursion
		sigma = getSigma(N, a_matrix, b_matrix, pi_matrix,
				sigma, numWords, obsIndex, a_ij_times_simga);
		chai = getChai(N, a_matrix, sigma, chai, numWords, a_ij_times_simga);
					System.out.println("chai: " + Arrays.deepToString(chai));
		//			System.out.println("numWords: " + numWords);
		//			System.out.println("sigma.length: " + sigma.length);
		//			System.out.println("sigma[0].length: " + sigma[0].length);

		// Step 3: termination step
		int q_star[]=new int[2];
		q_star = getQStar(sigma, chai, numWords, q_star);
		System.out.println("q_star[i]: The best sequence of states that will generate the given observation: " + Arrays.toString(q_star) );
		System.out.println("===============End of Viterabi Algorithm==============================");

		System.out.println("\n\n===============Start of Baum-Welch Algorithm==============================");


		for(int i=0; i<1; i++){
			runBaumWelch(N, M, a_matrix, b_matrix, pi_matrix, alpha, beta,
					numWords, obsIndex);
		}

		
		System.out.println("==============");

	}
	private static double[][] forwardProcedure(int N, double[][] a_matrix,
			double[][] b_matrix, double[] pi_matrix, int numWords,
			int[] obsIndex) {
		

		double[][] alpha;
		alpha = new double[numWords][N];
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
					alpha[t][l] += alpha[t-1][k] * a_matrix[k][l];
				}
			}
		}
		for(int t=1; t<numWords;t++){
			for(int j=0; j<N;j++){
				alpha[t][j]*= b_matrix[j][obsIndex[t]];
			}
		}
		System.out.println("alpha: "+ Arrays.deepToString(alpha));
		return alpha;
	}
	

	private static void runBaumWelch(int N, int M, double[][] a_matrix,
			double[][] b_matrix, double[] pi_matrix, double[][] alpha,
			double[][] beta, int numWords, int[] obsIndex) {
		double[][][] xi;
		double[][] gamma;
		xi = getXI(N, a_matrix, b_matrix, alpha, beta, numWords, obsIndex);
		System.out.println("xi: "+ Arrays.deepToString(xi));
		gamma = getGamma(N, numWords, alpha, beta);
		System.out.println("gamma: "+ Arrays.deepToString(gamma));

		for(int i=0; i<N;i++){
			pi_matrix[i] = gamma[0][i];
		}
		System.out.println("pi_matrix: " + Arrays.toString(pi_matrix));

		double[][] trained_a_matrix = new double[N][N];
		//			trained_a_matrix = a_matrix;
		for(int i=0; i<N;i++){
			for(int j=0;j<N;j++){
				for(int t=0; t<numWords-1; t++){
					trained_a_matrix[i][j] += xi[t][i][j];
				}
			}
		}

		for(int i=0;i<N;i++){
			for(int j=0; j<N;j++){
				double denom = 0.0;
				for(int t=0; t<numWords-1; t++){
					denom+=gamma[t][i];
				}
				if(denom==0.) denom = DENOMINATOR_CONSTANT;
				trained_a_matrix[i][j]/=denom;
				//					if(trained_a_matrix[i][j]>1.0) trained_a_matrix[i][j]=1.0;

			}
		}
		a_matrix = trained_a_matrix;
//		System.out.println("trained_a_matrix: " + Arrays.deepToString(trained_a_matrix));

		double[][] trained_b_matrix = new double[N][M];

		for(int i=0; i<N;i++){
			for(int m=0;m<M;m++){
				for(int t=0; t<numWords; t++){
					if(obsIndex[t]==m){
						trained_b_matrix[i][m]+=gamma[t][i];
					}
				}
			}
		}
		for(int i=0; i<N;i++){
			for(int m=0;m<M;m++){
				double denom =0.0;
				for(int t=0; t<numWords; t++){
					denom+= gamma[t][i];
				}
				if(denom==0.0)denom = DENOMINATOR_CONSTANT;
				trained_b_matrix[i][m]/=denom;
			}
		}
		b_matrix = trained_b_matrix;
//		System.out.println("trained_b_matrix: " + Arrays.deepToString(trained_b_matrix));
		alpha = getAlpha(N, trained_a_matrix, trained_b_matrix, pi_matrix, numWords, obsIndex);
		beta = getBeta(N, trained_a_matrix, trained_b_matrix, numWords, obsIndex);
//		System.out.println("trained_alpha_matrix: " + Arrays.deepToString(alpha));
//		System.out.println("trained_beta_matrix: " + Arrays.deepToString(beta));
		
		System.out.println("trained_a_matrix: " + Arrays.deepToString(a_matrix));
		System.out.println("trained_b_matrix: " + Arrays.deepToString(b_matrix));

	}

	private static double[][] getGamma(int N,  int numWords, double[][] alpha, double[][] beta) {
		double[][] gamma;
		gamma = new double[numWords][N];
		/*	
		//tutorial version
		for(int t=0; t<numWords;t++){
			for(int i=0; i<N;i++){
				for(int j=0;j<N;j++){
					gamma[t][i] += xi[t][i][j];
				}
			}
		}*/

		//wikipedia version
		for(int t=0; t<numWords;t++){
			for(int i=0; i<N;i++){
				gamma[t][i] = alpha[t][i]*beta[t][i];
				double sum=0.0;
				for(int j=0;j<N;j++){
					//System.out.println("j: " + j + ", i: " + i);
					sum+=alpha[t][j]*beta[t][j];
				}
				//System.out.println("sum: " +sum);
				if(sum==0.) sum = DENOMINATOR_CONSTANT;
				gamma[t][i]/=sum;
			}
		}
		return gamma;
	}

	private static double[][][] getXI(int N, double[][] a_matrix, double[][] b_matrix,
			double[][] alpha, double[][] beta, int numWords, int[] obsIndex) {
		double[][][] xi = new double[numWords][N][N];
		for(int t=0; t<numWords-1;t++){
			for(int i=0;i<N;i++){
				for(int j=0;j<N;j++){
					xi[t][i][j] = getXi_IJ_value(t, i, j , N, a_matrix, b_matrix, alpha, beta, xi, numWords, obsIndex);
				}
			}
		}
		return xi;
	}

	private static double getXi_IJ_value(int t, int state1, int state2, int N, double[][] a_matrix,
			double[][] b_matrix, double[][] alpha, double[][] beta,
			double[][][] xi, int numWords, int[] obsIndex) {
		xi[t][state1][state2] = alpha[t][state1] * a_matrix[state1][state2]*b_matrix[state2][obsIndex[t+1]] * beta[t+1][state2];


		//		System.out.println(" alpha[t][state1]  " +  alpha[t][state1] );
		//		System.out.println("b_matrix[state2][obsIndex[t+1]]: " + b_matrix[state2][obsIndex[t+1]]);
		//		System.out.println("beta[t+1][state2] " + beta[t+1][state2] +"\n");

		double denominator = 0.;
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				denominator+= alpha[t][i] * a_matrix[i][j] * b_matrix[j][obsIndex[t+1]] * beta[t+1][j];
			}
		}
/*
		//wiki version
		double denominator = 0.;
		for(int k=0;k<N;k++){
			denominator+= alpha[t][k] * beta[t][k];
		}
		//		System.out.println("xi[t][state1][state2]: " + xi[t][state1][state2]);
		//System.out.println("denominator: " + denominator);
*/
		if(denominator==0.) denominator = DENOMINATOR_CONSTANT;
		xi[t][state1][state2] /= denominator;
		//		System.out.println("xi[t][state1][state2]: " + xi[t][state1][state2]);
		return xi[t][state1][state2];
	}

	private static double[][] getBeta(int N, double[][] a_matrix,
			double[][] b_matrix, int numWords, int[] obsIndex) {
		double[][] beta;
		beta = new double[numWords][N];

		//initialization
		for(int i=0; i<N; i++){
			beta[numWords-1][i] = 1.0;
		}

		//		System.out.println("a_matrix : "+ Arrays.deepToString(a_matrix));
		//induction
		for(int t=numWords-2;t>=0;t--){
			for(int i=0; i<N;i++){
				double sum =0.0;
				for(int j=0;j<N;j++){
					sum+=a_matrix[i][j]*b_matrix[j][obsIndex[t+1]]*beta[t+1][j];
				}
				beta[t][i] = sum;
				//				System.out.println("beta: "+ Arrays.deepToString(beta));

			}
		}
		return beta;
	}

	private static int[] getQStar(double[][] sigma, double[][] chai,
			int numWords, int[] q_star) {
		double p_star = Arrays.stream(sigma[numWords-1]).max().getAsDouble();
		System.out.println("p_star: " + p_star);
		for(double e: sigma[numWords-1]){
			if(e==p_star) break;
			q_star[numWords-1]++;
		}

		//Step 4, Path.
		for(int i=2-2;i>=0;i--){
			q_star[i] =  (int) chai[i+1][(int)q_star[i+1]];
		}
		return q_star;
	}

	private static double[][] getChai(int N, double[][] a_matrix, double[][] sigma,
			double[][] chai, int numWords, double[] a_ij_times_simga) {
		//t=2 
		Random rand = new Random();
		for(int t=1;t<numWords;t++){
			for(int j=0;j<N;j++){
				for(int i=1;i<N;i++){
					a_ij_times_simga[i] = sigma[t-1][i]*a_matrix[i][j];
				}
				double highest_value = Arrays.stream(a_ij_times_simga).max().getAsDouble();
				/*
				if(highest_value!=0){
					int highest_index=0;
					for(double s:a_ij_times_simga){
						if(s==highest_value) break;
						highest_index++;
					}
					chai[t][j] = highest_index;
				} else {
					System.out.println("here");
					System.out.println("rand.nextInt(8): "+ rand.nextInt(8));
					chai[t][j] = (double)rand.nextInt(8);
				}
				*/
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

		//recursion step, t=2
		for(int t=1;t<numWords;t++){
			for(int j=0;j<N;j++){
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
					alpha[t][l] += alpha[t-1][k] * a_matrix[k][l];
				}
			}
		}
		for(int t=1; t<numWords;t++){
			for(int j=0; j<N;j++){
				alpha[t][j]*= b_matrix[j][obsIndex[t]];
			}
		}

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
