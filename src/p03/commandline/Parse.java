package p03.commandline;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class Parse {
	
	public static void main(String[] args) {
		
		int i = 0;
		int N = 0;
		int M = 0;
		int T = 0;
		
		double[][] matrixA = null;
		double[][] matrixB = null; 
		
		double[] matrixPI = null; 
		
		String currentLine = null; 
		String[] parsedLine = null;; 
		String[] vocabList = null; 
		
		BufferedReader buffer = null; 		
		
		if (args.length > 0) {
			if ((args[0].equals("./recognize") || args[0].equals("./statepath") || args[0].equals("./optimize")) && (args.length == 3 || args.length == 4)) {
				System.out.println("args[1]: " + args[1]); 
				System.out.println("args[2]: " + args[2]); 
				
				try {
					File input = new File(args[1]);
					File obs = new File(args[2]);	
					
					
					buffer = new BufferedReader(new FileReader(input));		
					
					while ((currentLine = buffer.readLine()) != null) {
						
						parsedLine = currentLine.split(" ");
						
						switch(i) {
						case 0: 
							N = Integer.parseInt(parsedLine[0]);
							M = Integer.parseInt(parsedLine[1]);
							T = Integer.parseInt(parsedLine[2]);
							
							matrixA = new double[N][N];
							matrixB = new double[N][M];
							matrixPI = new double[N];
							break; 
							
						case 1: 
							break; 
						case 2:
							vocabList = currentLine.split(" ");
							break;
						case 4:
							matrixA = extract_amatrix_Values(0, matrixA, currentLine.split(" "));
							break;
						case 5:
							matrixA = extract_amatrix_Values(1, matrixA, currentLine.split(" "));
							break;
						case 6:
							matrixA = extract_amatrix_Values(2, matrixA, currentLine.split(" "));
							break;
						case 7:
							matrixA = extract_amatrix_Values(3, matrixA, currentLine.split(" "));
							break;
						case 9:
							matrixB = extract_bmatrix_Values(0, matrixB, parsedLine);
							break;
						case 10:
							matrixB = extract_bmatrix_Values(1, matrixB, parsedLine);
							break;
						case 11:
							matrixB = extract_bmatrix_Values(2, matrixB, parsedLine);
							break;
						case 12:
							matrixB = extract_bmatrix_Values(3, matrixB, parsedLine);
							break;
						case 13: 
							break; 
						case 14:
							matrixPI[0] = Double.valueOf(parsedLine[0]);
							matrixPI[1] = Double.valueOf(parsedLine[1]);
							matrixPI[2] = Double.valueOf(parsedLine[2]);
							matrixPI[3] = Double.valueOf(parsedLine[3]);
							break;
						}

						i++;
					} // end of while
					
					// Initialization values so that the probability are non-zero
					for (int a = 0; a < 4; a++) {
						matrixA[a][0] = 0.23;
						matrixA[a][1] = 0.23;
						matrixA[a][2] = 0.21;
						matrixA[a][3] = 0.33;
						
						for (int b = 0; b < 4; b++) {
							matrixB[a][b] = 1 / 8; 
						}
						matrixPI[a] = 0.25; 
					}

					
					if (args[0].equals("./recognize")) {
						System.out.println("called recognize()");
//						recognize(obs,N,vocabList,matrixA,matrixB,matrixPI); 
					} else if (args[0].equals("./statepath")) { 
						
						System.out.println("called statepath()");
						find_statepath(); 
						
					} else if (args[0].equals("./optimize")) { 
						System.out.println("args[3]: " + args[3]);
						
						File output = new File(args[3]); 
						
						System.out.println("called optimize()");
						optimize(); 
						
					} else {
						System.out.println(args[0] +" does not match any existing methods.");
					}
				} catch (IOException error) { 
					error.printStackTrace(); 
				}
				
			} else { 
				System.out.println("Unrecognized command: " + args[0]);
			}
		}
	} 
	
	
	
	private static void optimize() {
		// TODO Auto-generated method stub
		
	}



	private static void find_statepath() {
		// TODO Auto-generated method stub
		
	}



//	private static void recognize(File obsFilename,int N, String[] list_of_vocabs,
//			double[][] a_matrix,double[][] b_matrix, double[] pi_matrix) throws NumberFormatException, IOException {
//		BufferedReader br=new BufferedReader(new FileReader(obsFilename));
//		String sCurrentLine =br.readLine() ;
//		int dataset_length = Integer.parseInt(sCurrentLine);
//
//		double[][] alpha;
//		double[][] beta;
//		double[][] sigma;
//		double[][] chai;
//		int numWords =0;
//
//		for(int dataSetIndex=0; dataSetIndex<dataset_length; dataSetIndex++){
//			if((sCurrentLine = br.readLine()) != null) numWords = Integer.parseInt(sCurrentLine);
//			String words[] =  br.readLine().split(" ");
//			int obsIndex[] = new int[numWords];
//			for(int j=0; j<numWords; j++){ 
//				obsIndex[j] = HMM_Algorithms.getIndex(words[j], list_of_vocabs);
//				//				System.out.println("obsIndex["+j+"]: " + obsIndex[j] );
//			}
//			alpha = HMM_Algorithms.getAlpha(N, a_matrix, b_matrix, pi_matrix, numWords,
//					obsIndex);
//			beta =  HMM_Algorithms.getBeta(N, a_matrix, b_matrix, numWords, obsIndex);
//			double Answer = 0.0;
//			for(int state=0; state<N; state++){ Answer += alpha[numWords-1][state]; }
////			System.out.println("alpha: "+ Arrays.deepToString(alpha));
////			System.out.println("beta: "+ Arrays.deepToString(beta));
//			System.out.println("The probablity of observing: "+Arrays.toString(words)+"=" + Answer);
////			
//		}
//		br.close();		
//	}



	private static double[][] extract_amatrix_Values(int i, double[][] a_matrix, String[] parseSpaceAvalues) {
		for (int k = 0; k < a_matrix[0].length; k++) {
			a_matrix[i][k] = Double.valueOf(parseSpaceAvalues[k]);
		}
		return a_matrix;
	}

	
	
	private static double[][] extract_bmatrix_Values(int i, double[][] b_matrix, String[] parseSpaceAvalues) {
		for (int k = 0; k < b_matrix[0].length; k++) {
			// System.out.println("Double.valueOf(parseSpaceAvalues[k]):  "
			// +Double.valueOf(parseSpaceAvalues[k]));
			b_matrix[i][k] = Double.valueOf(parseSpaceAvalues[k]);
		}
		return b_matrix;
	}
	
}