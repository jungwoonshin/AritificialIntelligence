package p03;

import java.util.Arrays;
import java.util.OptionalDouble;
import java.util.OptionalInt;

public class Testing {
	public static void main(String[] args) {
		int N = 2;
		double[][] a_matrix = {{0.69,0.3,0.01},{0.4,0.59,0.01}};
		double[][] b_matrix = {{0.5,0.4,0.1},{0.1,0.3,0.6}};
		int[] obs_index = {0,1,2};
		double[][] beta = getBeta( N, a_matrix,b_matrix, 3,obs_index);
		System.out.println("beta: " + Arrays.deepToString(beta));
/*
 * O = v1,v2,v3
 * PI = 0.6,0.4
 * A = [0.69,0.3,0.01;0.4,0.59,0.01]
 * B = [0.5,0.4,0.1;0.1,0.3,0.6]
 * ObsIndex = [0,1,2]
 */
	
//		states = ('Healthy', 'Fever')
//				end_state = 'E'
//				 
//				observations = ('normal', 'cold', 'dizzy')
//				 
//				start_probability = {'Healthy': 0.6, 'Fever': 0.4}
//				 
//				transition_probability = {
//				   'Healthy' : {'Healthy': 0.69, 'Fever': 0.3, 'E': 0.01},
//				   'Fever' : {'Healthy': 0.4, 'Fever': 0.59, 'E': 0.01},
//				   }
//				 
//				emission_probability = {
//				   'Healthy' : {'normal': 0.5, 'cold': 0.4, 'dizzy': 0.1},
//				   'Fever' : {'normal': 0.1, 'cold': 0.3, 'dizzy': 0.6},
//				   }
	
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
}
