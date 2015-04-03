package p03.commandline;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.OptionalDouble;

public class HMM_System {
    
    private static final double DENOMINATOR_CONSTANT = 10e-70;
    
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
        String[] stateList = null;
        
        BufferedReader buffer = null;
        
        if (args.length > 0) {
            if ((args[0].equals("./recognize") || args[0].equals("./statepath") || args[0].equals("./optimize")) && (args.length == 3 || args.length == 4)) {
//                System.out.println("args[1]: " + args[1]);
//                System.out.println("args[2]: " + args[2]);
                
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
                                stateList = new String[N];
                                break;
                                
                            case 1:
                                stateList = currentLine.split(" ");
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
                    }
                    
                    // Initialization values so that the probability are non-zero
                    for (int a = 0; a < 4; a++) {
                        matrixA[a][0] = 0.23;
                        matrixA[a][1] = 0.23;
                        matrixA[a][2] = 0.21;
                        matrixA[a][3] = 0.33;
                        
                        for (int b = 0; b < 4; b++) {
                            matrixB[a][b] = 1. / 8.;
                        }
                        matrixPI[a] = 0.25;
                    }
                    
                    
                    // Recognize Case
                    if (args[0].equals("./recognize")) {
                        if (args.length != 3) {
//                            System.out.println("calling recognize()\n");
                            
                            recognize(obs, N, vocabList, matrixA, matrixB, matrixPI);
                        }
                        
                        // Statepath Case
                    } else if (args[0].equals("./statepath")) {
                        statepath(obs,  N, vocabList, matrixA, matrixB, matrixPI,stateList);
                        
                        // Optimize Case
                    } else {
//                        System.out.println("args[3]: " + args[3]);
                        if(args[3]!=null) {File output = new File(args[3]);}
//                        System.out.println("calling optimize()\n");
                        
                        double[] probBefore = recognizeAux(obs, N, vocabList, matrixA, matrixB, matrixPI);
                        
                        ArrayList<double[][]> updatedMatrices = optimize(obs, vocabList, N, M, matrixA, matrixB, matrixPI);
                        
                        double[][] updatedMatrixA = updatedMatrices.get(0);
                        double[][] updatedMatrixB = updatedMatrices.get(1);
                        
                        double[] probAfter = recognizeAux(obs, N, vocabList, updatedMatrixA, updatedMatrixB, matrixPI);
                        
                        
                        for (int x = 0; x < probBefore.length; x++) {
                            System.out.println("Probability Before Optimization: " + probBefore[x]);
                            System.out.println("Probability After  Optimization: " + probAfter[x] + "\n");
                        }
                        
                    }
                    
                } catch (IOException error) {
                    error.printStackTrace();
                }
            } else {
                if (args.length != 3 && (args[0].equals("./recognize") || args[0].equals("./statepath"))) {
                    if (args.length < 3) {
                        System.out.println("Too few argument. recognize() or statepath() requires 3 args");
                    } else {
                        System.out.println("Too many argument. recognize() or statepath() requires 3 args");
                    }
                } else if (args.length != 4 && args[0].equals("./optimize")) {
                    if (args.length < 4) {
                        System.out.println("Too few arguments. optimize() requires 4 args");
                    } else {
                        System.out.println("Too many arguemnts. optimize() requires 4 args");
                    }
                } else {
                    System.out.println("Unrecognized command: " + args[0]);
                }
            }
        }
    }
    
    private static ArrayList<double[][]> optimize(File obsFilename, String[] vocablist, int N, int M, double[][] matrixA, double[][] matrixB, double[] matrixPI)
    throws NumberFormatException, IOException {
        BufferedReader br = new BufferedReader(new FileReader(obsFilename));
        String sCurrentLine = br.readLine();
        int dataset_length = Integer.parseInt(sCurrentLine);
        
        double[][] alpha;
        double[][] beta;
        double[][] sigma;
        double[][] chai;
        int numWords = 0;
        
        for (int dataSetIndex = 0; dataSetIndex < dataset_length; dataSetIndex++) {
            if ((sCurrentLine = br.readLine()) != null) {
                numWords = Integer.parseInt(sCurrentLine);
            }
            String words[] = br.readLine().split(" ");
            int obsIndex[] = new int[numWords];
            for (int j = 0; j < numWords; j++) {
                obsIndex[j] = getIndex(words[j], vocablist);
            }
            alpha = getAlpha(N, matrixA, matrixB, matrixPI, numWords, obsIndex);
            beta = getBeta(N, matrixA, matrixB, numWords, obsIndex);
            
            for (int i = 0; i < 1; i++) {
                runBaumWelch(N, M, matrixA, matrixB, matrixPI, alpha, beta, numWords, obsIndex);
            }
        }
//        System.out.println("Updated A: " + Arrays.deepToString(matrixA));
//        System.out.println("Updated B: " + Arrays.deepToString(matrixB));
        
        
        ArrayList<double[][]> updatedMatrices = new ArrayList<double[][]>();
        updatedMatrices.add(matrixA);
        updatedMatrices.add(matrixB);
        
        br.close();
        
        return updatedMatrices;
 
        // Run the forward algorithm again with the updated matrices
        // Output the before and after optimization probabilities
    }
    
    private static void statepath(File obsFilename, int N, String[] list_of_vocabs, double[][] a_matrix, double[][] b_matrix, double[] pi_matrix, String[] stateList) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(obsFilename));
        String sCurrentLine = br.readLine();
        int dataset_length = Integer.parseInt(sCurrentLine);
        
        double[][] alpha;
        double[][] beta;
        double[][] sigma;
        double[][] chai;
        int numWords = 0;
        
        for (int dataSetIndex = 0; dataSetIndex < dataset_length; dataSetIndex++) {
            if ((sCurrentLine = br.readLine()) != null) {
                numWords = Integer.parseInt(sCurrentLine);
            }
            String words[] = br.readLine().split(" ");
            int obsIndex[] = new int[numWords];
            for (int j = 0; j < numWords; j++) {
                obsIndex[j] = getIndex(words[j], list_of_vocabs);
            }
            alpha = getAlpha(N, a_matrix, b_matrix, pi_matrix, numWords, obsIndex);
            beta = getBeta(N, a_matrix, b_matrix, numWords, obsIndex);
            runViterabi(N, list_of_vocabs, a_matrix, b_matrix, pi_matrix, numWords, obsIndex, stateList);
        }
        br.close();
    }
    
    private static void recognize(File obsFilename, int N, String[] list_of_vocabs, double[][] a_matrix, double[][] b_matrix, double[] pi_matrix)
    throws NumberFormatException, IOException {
        BufferedReader br = new BufferedReader(new FileReader(obsFilename));
        String sCurrentLine = br.readLine();
        int dataset_length = Integer.parseInt(sCurrentLine);
        
        double[][] alpha;
        double[][] beta;
        double[][] sigma;
        double[][] chai;
        int numWords = 0;
        
        for (int dataSetIndex = 0; dataSetIndex < dataset_length; dataSetIndex++) {
            if ((sCurrentLine = br.readLine()) != null) {
                numWords = Integer.parseInt(sCurrentLine);
            }
            String words[] = br.readLine().split(" ");
            int obsIndex[] = new int[numWords];
            for (int j = 0; j < numWords; j++) {
                obsIndex[j] = getIndex(words[j], list_of_vocabs);
            }
            alpha = getAlpha(N, a_matrix, b_matrix, pi_matrix, numWords, obsIndex);
            beta = getBeta(N, a_matrix, b_matrix, numWords, obsIndex);
            double Answer = 0.0;
            for (int state = 0; state < N; state++) {
                Answer += alpha[numWords - 1][state];
            }
            
            System.out.println(Answer + " " + Arrays.toString(words));
        }
        
        br.close();
    }
    
    private static double[] recognizeAux(File obsFilename, int N, String[] list_of_vocabs, double[][] a_matrix, double[][] b_matrix, double[] pi_matrix)
    throws NumberFormatException, IOException {
        BufferedReader br = new BufferedReader(new FileReader(obsFilename));
        String sCurrentLine = br.readLine();
        int dataset_length = Integer.parseInt(sCurrentLine);
        
        double[][] alpha;
        double[][] beta;
        double[][] sigma;
        double[][] chai;
        int numWords = 0;
        
        double[] ret = new double[dataset_length];
        
        for (int dataSetIndex = 0; dataSetIndex < dataset_length; dataSetIndex++) {
            if ((sCurrentLine = br.readLine()) != null) {
                numWords = Integer.parseInt(sCurrentLine);
            }
            String words[] = br.readLine().split(" ");
            int obsIndex[] = new int[numWords];
            for (int j = 0; j < numWords; j++) {
                obsIndex[j] = getIndex(words[j], list_of_vocabs);
            }
            alpha = getAlpha(N, a_matrix, b_matrix, pi_matrix, numWords, obsIndex);
            beta = getBeta(N, a_matrix, b_matrix, numWords, obsIndex);
            
            double Answer = 0.0;
            for (int state = 0; state < N; state++) {
                Answer += alpha[numWords - 1][state];
            }
            
            ret[dataSetIndex] = Answer;
            
        }
        
        br.close();
        return ret;
    }
    
    
    public static void runViterabi(int N, String[] list_of_vocabs, double[][] a_matrix, double[][] b_matrix, double[] pi_matrix, int numWords, int[] obsIndex, String[] stateList) {
        double[][] sigma;
        double[][] chai;
        sigma = new double[numWords][N];
        chai = new double[numWords][N];
        double[] a_ij_times_sigma = new double[N];
        
        // Step 1 & 2: Recursion
        sigma = getSigma(N, a_matrix, b_matrix, pi_matrix, sigma, numWords, obsIndex, a_ij_times_sigma);
        chai = getChai(N, a_matrix, sigma, chai, numWords, a_ij_times_sigma);
        
        // Step 3: termination step
        int q_star[] = new int[numWords];
        q_star = getQStar(sigma, chai, numWords, q_star);
        
        
        System.out.print("Generated the Best Plausible State Sequence: ");
        
        for (int i = 0; i < q_star.length; i++) {
            System.out.print(stateList[q_star[i]] + " ");
        }
        System.out.println("\n");
    }
    
    public static void runBaumWelch(int N, int M, double[][] a_matrix, double[][] b_matrix, double[] pi_matrix, double[][] alpha, double[][] beta,
                                    int numWords, int[] obsIndex) {
        double[][][] xi;
        double[][] gamma;
        xi = getXI(N, a_matrix, b_matrix, alpha, beta, numWords, obsIndex);
        gamma = getGamma(N, numWords, alpha, beta);
        
        for (int i = 0; i < N; i++) {
            pi_matrix[i] = gamma[0][i];
        }
        
        double[][] trained_a_matrix = new double[N][N];
        
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                for (int t = 0; t < numWords - 1; t++) {
                    trained_a_matrix[i][j] += xi[t][i][j];
                }
            }
        }
        
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                double denom = 0.0;
                for (int t = 0; t < numWords - 1; t++) {
                    denom += gamma[t][i];
                }
                if (denom == 0.) {
                    denom = DENOMINATOR_CONSTANT;
                }
                trained_a_matrix[i][j] /= denom;
                
            }
        }
        a_matrix = trained_a_matrix;
        
        double[][] trained_b_matrix = new double[N][M];
        
        for (int i = 0; i < N; i++) {
            for (int m = 0; m < M; m++) {
                for (int t = 0; t < numWords; t++) {
                    if (obsIndex[t] == m) {
                        trained_b_matrix[i][m] += gamma[t][i];
                    }
                }
            }
        }
        for (int i = 0; i < N; i++) {
            for (int m = 0; m < M; m++) {
                double denom = 0.0;
                for (int t = 0; t < numWords; t++) {
                    denom += gamma[t][i];
                }
                if (denom == 0.0) {
                    denom = DENOMINATOR_CONSTANT;
                }
                trained_b_matrix[i][m] /= denom;
            }
        }
        b_matrix = trained_b_matrix;
        alpha = getAlpha(N, trained_a_matrix, trained_b_matrix, pi_matrix, numWords, obsIndex);
        beta = getBeta(N, trained_a_matrix, trained_b_matrix, numWords, obsIndex);
    }
    
    public static double[][] getGamma(int N, int numWords, double[][] alpha, double[][] beta) {
        double[][] gamma;
        gamma = new double[numWords][N];
        
        // Wikipedia version
        for (int t = 0; t < numWords; t++) {
            for (int i = 0; i < N; i++) {
                gamma[t][i] = alpha[t][i] * beta[t][i];
                double sum = 0.0;
                for (int j = 0; j < N; j++) {
                    sum += alpha[t][j] * beta[t][j];
                }
                if (sum == 0.) {
                    sum = DENOMINATOR_CONSTANT;
                }
                gamma[t][i] /= sum;
            }
        }
        return gamma;
    }
    
    public static double[][][] getXI(int N, double[][] a_matrix, double[][] b_matrix, double[][] alpha, double[][] beta, int numWords, int[] obsIndex) {
        double[][][] xi = new double[numWords][N][N];
        for (int t = 0; t < numWords - 1; t++) {
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    xi[t][i][j] = getXi_IJ_value(t, i, j, N, a_matrix, b_matrix, alpha, beta, xi, numWords, obsIndex);
                }
            }
        }
        return xi;
    }
    
    public static double getXi_IJ_value(int t, int state1, int state2, int N, double[][] a_matrix, double[][] b_matrix, double[][] alpha,
                                        double[][] beta, double[][][] xi, int numWords, int[] obsIndex) {
        xi[t][state1][state2] = alpha[t][state1] * a_matrix[state1][state2] * b_matrix[state2][obsIndex[t + 1]] * beta[t + 1][state2];
        
        // Wiki version
        double denominator = 0.;
        for (int k = 0; k < N; k++) {
            denominator += alpha[t][k] * beta[t][k];
        }
        
        if (denominator == 0.) {
            denominator = DENOMINATOR_CONSTANT;
        }
        xi[t][state1][state2] /= denominator;
        
        return xi[t][state1][state2];
    }
    
    public static double[][] getBeta(int N, double[][] a_matrix, double[][] b_matrix, int numWords, int[] obsIndex) {
        double[][] beta;
        beta = new double[numWords][N];
        
        // Initialization
        for (int i = 0; i < N; i++) {
            beta[numWords - 1][i] = 1.0;
        }
        
        // Induction
        for (int t = numWords - 2; t >= 0; t--) {
            for (int i = 0; i < N; i++) {
                double sum = 0.0;
                for (int j = 0; j < N; j++) {
                    sum += a_matrix[i][j] * b_matrix[j][obsIndex[t + 1]] * beta[t + 1][j];
                }
                beta[t][i] = sum;
            }
        }
        return beta;
    }
    
    public static int[] getQStar(double[][] sigma, double[][] chai, int numWords, int[] q_star) {
        double p_star = Arrays.stream(sigma[numWords - 1]).max().getAsDouble();
        
        for (double e : sigma[numWords - 1]) {
            if (e == p_star) {
                break;
            }
            q_star[numWords - 1]++;
        }
        
        // Step 4:  Path
        for (int i = numWords - 2; i >= 0; i--) {
            q_star[i] = (int) chai[i + 1][(int) q_star[i + 1]];
        }
        return q_star;
    }
    
    public static double[][] getChai(int N, double[][] a_matrix, double[][] sigma, double[][] chai, int numWords, double[] a_ij_times_simga) {
        
        for (int t = 1; t < numWords; t++) {
            for (int j = 0; j < N; j++) {
                for (int i = 1; i < N; i++) {
                    a_ij_times_simga[i] = sigma[t - 1][i] * a_matrix[i][j];
                }
                double highest_value = Arrays.stream(a_ij_times_simga).max().getAsDouble();
                
                int highest_index = 0;
                for (double s : a_ij_times_simga) {
                    if (s == highest_value) {
                        break;
                    }
                    highest_index++;
                }
                chai[t][j] = highest_index;
            }
        }
        return chai;
    }
    
    public static double[][] getSigma(int N, double[][] a_matrix, double[][] b_matrix, double[] pi_matrix, double[][] sigma, int numWords,
                                      int[] obsIndex, double[] a_ij_times_simga) {
        // Initialization Step: t = 1
        for (int state = 0; state < N; state++) {
            sigma[0][state] = pi_matrix[state] * b_matrix[state][obsIndex[0]];
        }
        
        // Recursion Step, t = 2
        for (int t = 1; t < numWords; t++) {
            for (int j = 0; j < N; j++) {
                for (int state = 0; state < N; state++) {
                    a_ij_times_simga[state] = sigma[t - 1][state] * a_matrix[state][j];
                }
                OptionalDouble highest = Arrays.stream(a_ij_times_simga).max();
                sigma[t][j] = highest.getAsDouble() * b_matrix[j][obsIndex[t]];
            }
        }
        return sigma;
    }
    
    /*
     * Executes Forward Procedure and returns the alpha array. Size of alpha is
     * T*N where T is the number of observed time sequence and N is the number
     * of possible states.
     */
    public static double[][] getAlpha(int N, double[][] a_matrix, double[][] b_matrix, double[] pi_matrix, int numWords, int[] obsIndex) {
        double[][] alpha = new double[numWords][N];
        // Step : t = 1
        for (int state = 0; state < N; state++) {
            alpha[0][state] = pi_matrix[state] * b_matrix[state][obsIndex[0]];
        }
        
        // Step : t = 2
        for (int t = 1; t < numWords; t++) {
            for (int l = 0; l < N; l++) {
                for (int k = 0; k < N; k++) {
                    alpha[t][l] += alpha[t - 1][k] * a_matrix[k][l];
                }
            }
        }
        
        for (int t = 1; t < numWords; t++) {
            for (int j = 0; j < N; j++) {
                alpha[t][j] *= b_matrix[j][obsIndex[t]];
            }
        }
        
        return alpha;
    }
    
    /*
     * From each observed token sequence, finds the index that matches the
     * vocabulary list. It can return [0,m] where m is the number of possible
     * outputs.
     */
    public static int getIndex(String word, String vocab[]) {
        int i = 0;
        while (i < vocab.length && !vocab[i].equalsIgnoreCase(word)) {
            i++;
        }
        return i;
    }
    
    /*
     * Reads matrix A and returns N*N array. N is the number of states
     */
    public static double[][] extract_amatrix_Values(int i, double[][] a_matrix, String[] parseSpaceAvalues) {
        for (int k = 0; k < a_matrix[0].length; k++) {
            a_matrix[i][k] = Double.valueOf(parseSpaceAvalues[k]);
        }
        return a_matrix;
    }
    
    /*
     * Read matrix B and returns N*M array. N is number of states. M is the
     * number of possible outputs.
     */
    public static double[][] extract_bmatrix_Values(int i, double[][] b_matrix, String[] parseSpaceAvalues) {
        for (int k = 0; k < b_matrix[0].length; k++) {
            b_matrix[i][k] = Double.valueOf(parseSpaceAvalues[k]);
        }
        return b_matrix;
    }
    
}