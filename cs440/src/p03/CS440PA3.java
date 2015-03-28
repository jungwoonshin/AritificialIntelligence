package p03;

/*
 * Hidden Markov Models and Natural Language Processing
 * CS440 Programming Assignment 3
 * 
 * Code written by:
 * Adam Engel
 * Chris Hall
 * Tim Duffy
 * 
 */

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

public class CS440PA3 {

	public static void main(String args[]) throws FileNotFoundException /*throws FileNotFoundException */{
		String line1 = "/Users/jungwoonshin/git/cs440/cs440/src/p03/sentence.hmm";
		String line2 = "/Users/jungwoonshin/git/cs440/cs440/src/p03/example1.obs";
		String line3 = "/Users/jungwoonshin/git/cs440/cs440/src/p03/example2.obs";
		recognize(line1,line2 );
		recognize(line1,line3 );
//		System.out.println("Please enter a command with the appropriate arguments:");
//		System.out.print(System.getProperty("user.name")+"$ ");
//
//		Scanner console = new Scanner(System.in);
//		String line[] = console.nextLine().split(" ");
//
//		while(!line[0].equalsIgnoreCase("q") && !line[0].equalsIgnoreCase("quit")){
//
//			try
//			{
//
//				// If the user selects the 'recognize' function (forward algorithm)
//				if(line[0].equalsIgnoreCase("./recognize") || line[0].equalsIgnoreCase("recognize")){
//					if(line.length != 3){
//						System.out.println("Incorrect number of arguments. \nUsage is ./recognize <HMM> <observation sets>");
//					}
//					else{
//						recognize(line[1], line[2]);
//						System.out.println();
//					}
//				}
//				else if(line[0].equalsIgnoreCase("./statepath") || line[0].equalsIgnoreCase("statepath")){
//					if(line.length != 3){
//						System.out.println("Incorrect number of arguments. \nUsage is ./statepath <HMM> <observation sets>");
//					}
//					else{
//						statepath(line[1], line[2]);
//					}
//				}
//				else if(line[0].equalsIgnoreCase("./optimize") || line[0].equalsIgnoreCase("optimize")){
//					if(line.length != 4){
//						System.out.println("Incorrect number of arguments. \nUsage is ./optimize <HMM> <observation sets> <output file>");
//					}
//					else{
//						recognize(line[1], line[2]);			// First, print output for original HMM
//						optimize(line[1], line[2], line[3]);	// Apply the Baum-Welch algorithm 1 time
//					}
//				}
//				else{
//					System.out.println("\nIncorrect command. Usage is:");
//					System.out.println("./recognize <HMM> <observation sets>");
//					System.out.println("./statepath <HMM> <observation sets>");
//					System.out.println("./optimize <HMM> <observation sets> <output file>");
//					System.out.println();
//				}
//			}
//			catch(IOException e){
//				System.out.println("One of your files was not found. Check the paths and try again.");
//			}
//
//			// Print out username like a real command prompt
//			System.out.print(System.getProperty("user.name")+"$ ");
//			line = console.nextLine().split(" ");
//		}  
//
//		System.out.println("+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+");
	}

	/*
	 * The recognize function takes a file containing information about the Hidden Markov Model and a 
	 * file containing observation sets as its parameters. It then uses the 'forward' part of the
	 * 'forward-backward' algorithm to determine the probability of seeing a particular sequence of 
	 * observation symbols, and returns an array of these probabilities for each set in the file
	 * containing observation sets.
	 */
	public static void recognize(String HMMFile, String obSetsFile) throws FileNotFoundException{

		int numStates;
		int numSymbols;
		int T;
		String vocab[];
		State HMM[];

		FileReader f = new FileReader(HMMFile);
		Scanner hmmIn = null;

		// Read information from HMM file
		try {

			hmmIn = new Scanner(new BufferedReader(f));
			hmmIn.useLocale(Locale.US);

			// The first line contains integers N (number of states), M (number of observation symbols)
			// and T (number of time steps or length of observation sequences).
			String str[] = hmmIn.nextLine().split(" ");
			numStates = Integer.parseInt(str[0]);
			numSymbols = Integer.parseInt(str[1]);
			T = Integer.parseInt(str[2]);

			// The second line contains the name of the states
			String stateNames[] = hmmIn.nextLine().split(" ");

			// The third line contains the vocabulary to be used in the observation sentences
			vocab = new String[numSymbols]; 
			vocab = hmmIn.nextLine().split(" ");

			hmmIn.nextLine();							// skip over the text "a:"
			String A[] = new String[numStates];		// A holds the transition probabilities to each state
			for(int i=0; i<numStates; i++){ A[i] = hmmIn.nextLine(); }

			hmmIn.nextLine();							// skip over the text "b:"
			String B[] = new String[numStates];		// B holds the observation probabilities of vocab words in each state
			for(int i=0; i<numStates; i++){ B[i] = hmmIn.nextLine(); }

			hmmIn.nextLine();							// skip over the text "pi:"
			String pi[] = hmmIn.nextLine().split(" ");				// Starting probabilities

			// Initialize each state 
			HMM = new State[numStates];
			for(int i=0; i<numStates; i++){
				HMM[i] = new State(stateNames[i], i, A[i], B[i], Double.parseDouble(pi[i]));
			}
		}finally { hmmIn.close(); }

		// Read in observation sets
		FileReader f2 = new FileReader(obSetsFile);
		Scanner obsIn = null;

		try{

			obsIn = new Scanner(new BufferedReader(f2));
			obsIn.useLocale(Locale.US);

			int numSets = Integer.parseInt(obsIn.nextLine());

			// The Forward algorithm
			for(int obsSet=0; obsSet<numSets; obsSet++){

//				if(obsSet > 0) { System.out.println(); }

				int numWords = Integer.parseInt(obsIn.nextLine());
				String words[] = obsIn.nextLine().split(" ");
				int seqLen;

				if(numWords >= T){ seqLen = numWords+1; }
				else { seqLen = T; }

				// Initialization and Induction Steps are done by the forwardAlgo()
				double probs[][] = forwardAlgo(HMM, vocab, words, numWords, seqLen);

				// Termination Step
				double seqProb = 0.0;
				for(int state=0; state<numStates; state++){ seqProb += probs[seqLen-1][state]; }
				
				double seqProb1 = 0.0;
				for(int state=0; state<numStates; state++){ seqProb1 += probs[seqLen][state]; }
				
				double seqProb2 = 0.0;
				for(int state=0; state<numStates; state++){ seqProb2 += probs[seqLen-2][state]; }

				if(seqProb > 0.0){	// If the probability is non-zero, format it nicely
					DecimalFormat df = new DecimalFormat("#.######");
//					System.out.print(df.format(seqProb));
				}
				else{	// Otherwise print as 0.0
//					System.out.print(seqProb);
				}
				double probs2[][] = backwardAlgo(HMM, vocab, words, numWords, seqLen);
//				System.out.println("alpha: " + Arrays.deepToString(probs));
				System.out.println("beta: " + Arrays.deepToString(probs2));
			}
		}finally { obsIn.close(); }
	}

	/*
	 * The statepath function takes a file containing information about the Hidden Markov Model
	 * and a file containing observation sets as its parameters. It then uses the Viterbi algorithm
	 * to determine (and output) the path of states with the highest probability 
	 */
	public static void statepath(String HMMFile, String obSetsFile) throws FileNotFoundException{

		int numStates;
		int numSymbols;
		int T;
		String vocab[];
		State HMM[];

		FileReader f = new FileReader(HMMFile);
		Scanner hmmIn = null;

		// Read information from HMM file
		try {

			hmmIn = new Scanner(new BufferedReader(f));
			hmmIn.useLocale(Locale.US);

			// The first line contains integers N (number of states), M (number of observation symbols)
			// and T (number of time steps or length of observation sequences).
			String str[] = hmmIn.nextLine().split(" ");
			numStates = Integer.parseInt(str[0]);
			numSymbols = Integer.parseInt(str[1]);
			T = Integer.parseInt(str[2]);

			// The second line contains the name of the states
			String stateNames[] = hmmIn.nextLine().split(" ");

			// The third line contains the vocabulary to be used in the observation sentences
			vocab = new String[numSymbols]; 
			vocab = hmmIn.nextLine().split(" ");

			hmmIn.nextLine();						// skip over the text "a:"
			String A[] = new String[numStates];		// A holds the transition probabilities to each state
			for(int i=0; i<numStates; i++){ A[i] = hmmIn.nextLine(); }

			hmmIn.nextLine();						// skip over the text "b:"
			String B[] = new String[numStates];		// B holds the observation probabilities of vocab words in each state
			for(int i=0; i<numStates; i++){ B[i] = hmmIn.nextLine(); }

			hmmIn.nextLine();							// skip over the text "pi:"
			String pi[] = hmmIn.nextLine().split(" ");	// Starting probabilities

			// Initialize each state 
			HMM = new State[numStates];
			for(int i=0; i<numStates; i++){
				HMM[i] = new State(stateNames[i], i, A[i], B[i], Double.parseDouble(pi[i]));
			}
		}finally { hmmIn.close(); }

		// Read in observation sets
		FileReader f2 = new FileReader(obSetsFile);
		Scanner obsIn = null;

		try{

			obsIn = new Scanner(new BufferedReader(f2));
			obsIn.useLocale(Locale.US);

			int numSets = Integer.parseInt(obsIn.nextLine());

			for(int obsSet=0; obsSet<numSets; obsSet++){

				int numWords = Integer.parseInt(obsIn.nextLine());
				String words[] = obsIn.nextLine().split(" ");
				int seqLen;

				if(numWords >= T){ seqLen = numWords+1; }
				else { seqLen = T; }

				// Determine the index of each observation symbol (so we can find it in the A and B matrices)
				int obsIndex[] = new int[numWords+1];
				for(int i=0; i<numWords; i++){ obsIndex[i+1] = getIndex(words[i], vocab); }

				// Veterbi v holds the arrays p(probabilities along most likely path) and q(most likely path of states)
				Viterbi v = viterbiAlgo(HMM, vocab, words, numWords, seqLen);

				// Print probability of most likely path
				if(v.p[seqLen-1] > 0.0){		// If the probability is non-zero, format it nicely
					DecimalFormat df = new DecimalFormat("#.######");
					System.out.print(df.format(v.p[seqLen-1]) + " ");
				}
				else{					// Otherwise print as 0.0
					System.out.print(v.p[seqLen-1] + " ");
				}

				// Read out path
				if(v.p[seqLen -1] > 0.0){
					for(int t=1; t<seqLen; t++){
						System.out.print( HMM[ v.q[t] ].name + " ");
					}
				}
				System.out.println();
			}
		}finally { obsIn.close(); }
	}

	public static void optimize(String HMMFile, String obSetsFile, String outputFile) throws FileNotFoundException{

		int numStates;
		int numSymbols;
		int T;
		String vocab[];
		State HMM[];
		DecimalFormat df = new DecimalFormat("#.######");
		double newPi[] = null;
		double newA[][] = null;
		double newB[][] = null;

		FileReader f = new FileReader(HMMFile);
		Scanner hmmIn = null;

		// Read information from HMM file
		try {

			hmmIn = new Scanner(new BufferedReader(f));
			hmmIn.useLocale(Locale.US);

			// The first line contains integers N (number of states), M (number of observation symbols)
			// and T (number of time steps or length of observation sequences).
			String str[] = hmmIn.nextLine().split(" ");
			numStates = Integer.parseInt(str[0]);
			numSymbols = Integer.parseInt(str[1]);
			T = Integer.parseInt(str[2]);

			// The second line contains the name of the states
			String stateNames[] = hmmIn.nextLine().split(" ");

			// The third line contains the vocabulary to be used in the observation sentences
			vocab = new String[numSymbols];
			vocab = hmmIn.nextLine().split(" ");

			hmmIn.nextLine();							// skip over the text "a:"
			String A[] = new String[numStates];		// A holds the transition probabilities to each state
			for(int i=0; i<numStates; i++){ A[i] = hmmIn.nextLine(); }

			hmmIn.nextLine();							// skip over the text "b:"
			String B[] = new String[numStates];		// B holds the observation probabilities of vocab words in each state
			for(int i=0; i<numStates; i++){ B[i] = hmmIn.nextLine(); }

			hmmIn.nextLine();							// skip over the text "pi:"
			String pi[] = hmmIn.nextLine().split(" ");				// Starting probabilities

			// Initialize each state 
			HMM = new State[numStates];
			for(int i=0; i<numStates; i++){
				HMM[i] = new State(stateNames[i], i, A[i], B[i], Double.parseDouble(pi[i]));
			}
		}finally { hmmIn.close(); }

		// Read in observation sets
		FileReader f2 = new FileReader(obSetsFile);
		Scanner obsIn = null;

		try{

			obsIn = new Scanner(new BufferedReader(f2));
			obsIn.useLocale(Locale.US);

			int numSets = Integer.parseInt(obsIn.nextLine());

			// The Baum-Welsh Algorithm
			for(int obsSet=0; obsSet<numSets; obsSet++){

				int seqLen;
				int numWords = Integer.parseInt(obsIn.nextLine());
				String words[] = obsIn.nextLine().split(" ");

				if(numWords >= T){ seqLen = numWords+1; }
				else { seqLen = T; }

				double forwardProbs[][] = forwardAlgo(HMM, vocab, words, numWords, seqLen);
				double backwardProbs[][] = backwardAlgo(HMM, vocab, words, numWords, seqLen);
				
				/*
				System.out.println("alpha:");
				for(int i = 0; i<numStates; i++){
					for(int t = 1; t<=seqLen; t++){
						System.out.print(forwardProbs[t][i] + "\t");
					}System.out.println();
				}
				
				System.out.println("beta:");
				for(int i = 0; i<numStates; i++){
					for(int t = 1; t<=seqLen; t++){
						System.out.print(backwardProbs[t][i] + "\t");
					}System.out.println();
				} */

				// Determine the index of each observation symbol (so we can find it in the A and B matrices)
				int obsIndex[] = new int[seqLen+1];
				for(int i=0; i<numWords; i++)
				{ 
					obsIndex[i+1] = getIndex(words[i], vocab);
				}

				// xi[t][i][j] holds the probability of being in state <i> at time <t> and transitioning to state j at time <t+1>
				double xi[][][] = new double[seqLen+1][numStates][numStates];
				
				// Calculate xi tables
				for(int t=0; t<seqLen; t++){
					
					double seqProb = 0.0;
					for(int i=0; i<numStates; i++){
						seqProb += forwardProbs[t][i] * backwardProbs[t][i];
					}
					
					for(int i=0; i<numStates; i++){
						for(int j=0; j<numStates; j++){
							xi[t][i][j]= (forwardProbs[t][i] * HMM[i].A[j] * HMM[j].B[obsIndex[t+1]] * backwardProbs[t+1][j]) / seqProb;
						}
					}
				}

				// gamma[t][i] holds the probability of being in state <i> at time <t>
				double gamma[][] = new double[seqLen+1][numStates];

				// Calculate gamma tables
				for(int t=1; t<seqLen; t++){
					
					double seqProb = 0.0;
					for(int i=0; i<numStates; i++){
						seqProb += forwardProbs[t][i] * backwardProbs[t][i];
					}
					
					for(int state=0; state<numStates; state++){
						gamma[t][state] = (forwardProbs[t][state] * backwardProbs[t][state]) / seqProb;						
					}
				}
				
				/*
				System.out.println("gammas:");
				for(int i = 0; i<numStates; i++){
					for(int t = 1; t<seqLen; t++){
						System.out.print(gamma[t][i] + "\t");
					}System.out.println();
				}
				*/
				
				// Re-estimation formulas
				newPi = new double[numStates];
				newA = new double[numStates][numStates];
				newB = new double[numStates][numSymbols];

				// Calculate new pi
				for(int i=0; i<numStates; i++){
					newPi[i] = gamma[1][i];
				}

				// Calculate new transition probabilities
				for(int i=0; i<numStates; i++){
					for(int j=0; j<numStates; j++){

						double num = 0.0;
						double denominator = 0.0;
						
						for(int t=1; t<seqLen; t++){
							num += xi[t][i][j];
							denominator += gamma[t][i];
						}
						
						if (denominator > 0) { newA[i][j] = num / denominator; }
						else { newA[i][j] = HMM[i].A[j]; }
					}
				}

				// Calculate new observation probabilities
				for(int i=0; i<numStates; i++){
					for(int k=0; k<vocab.length; k++){

						double num = 0.0;
						double denom = 0.0;

						for(int t=1; t<seqLen; t++){
							
							denom += gamma[t][i];
							
							// only add to numerator if o<t> = k
							if(obsIndex[t] == k) { num += gamma[t][i]; }
						}

						if (denom > 0.0){ 
							newB[i][k] = num/denom;
						}
						else
						{ newB[i][k] = HMM[i].B[k]; }
					}
				}

				// Write new estimates to file
				try{
					// Create file 
					FileWriter fstream = new FileWriter(outputFile);
					BufferedWriter out = new BufferedWriter(fstream);

					// The first line contains integers N (number of states), M (number of observation symbols)
					// and T (number of time steps or length of observation sequences).
					out.write(numStates + " " + numSymbols + " " + words.length + "\n");

					// The second line contains the name of the states
					for(int i=0; i<numStates; i++){ out.write(HMM[i].name + " "); }
					out.write("\n");

					// The third line contains the vocabulary to be used in the observation sentences
					for(int i=0; i<numSymbols; i++){ out.write(vocab[i] + " "); }
					out.write("\n");

					// Write the transition probability matrix A
					out.write("a:\n");
					for(int i=0; i<numStates; i++){
						for(int j=0; j<numStates; j++){ 
							out.write( df.format(newA[i][j]) + " ");
						}
						out.write("\n");
					}

					// Write the observation probability matrix B
					out.write("b:\n");
					for(int i=0; i<numStates; i++){
						for(int j=0; j<numSymbols; j++){ out.write(df.format(newB[i][j]) + " "); }
						out.write("\n");
					}

					// Write the starting probabilities of each state
					out.write("pi:\n");
					for(int i=0; i<numStates; i++){ out.write(df.format(newPi[i]) + " "); }
					out.write("\n");

					//Close the output stream
					out.close();
				}catch (Exception e){//Catch exception if any
					System.err.println("Error: " + e.getMessage());
				}
			}
		}finally { obsIn.close(); }
		
		// Finally, output the probability of the sequence, using the optimized HMM
		System.out.print("  ");
		recognize(outputFile, obSetsFile);
		System.out.println();
	}

	/*
	 * Given an array of hidden Markov states, an array of vocabulary words, an observation set, the number of
	 * words in that set, and the length of an observation sequence, this algorithm calculates the forward
	 * probability of observing the given sequence. It then returns the table of forward probabilities for
	 * each state and time t.   
	 */
	public static double[][] forwardAlgo(State HMM[], String vocab[], String words[], int numWords, int seqLen){

		int numStates = HMM.length;

		// Determine the index of each observation symbol (so we can find it in the A and B matrices)
		int obsIndex[] = new int[seqLen+1];
		for(int i=0; i<numWords; i++){ obsIndex[i+1] = getIndex(words[i], vocab); }

		// probs[t][i] holds the probability of an observation sequence of length <t> ending in state <i>
		double probs[][] = new double[seqLen+1][numStates];

		// Initialize pi probabilities
		for(int state=0; state<numStates; state++){
			probs[0][state] = HMM[state].pi;
		}

		// Initialization Step: t = 1
		for(int state=0; state<numStates; state++){
			probs[1][state] =  HMM[state].pi * HMM[state].B[obsIndex[1]];
		}

		// Induction Step: t = 2.. seqLen-1
		for(int t=1; t<seqLen; t++){
			for(int j=0; j<numStates; j++){

				double sum = 0.0;
				for(int i=0; i<numStates; i++){
					sum += probs[t][i] * HMM[i].A[j];
				}
				
				if(HMM[j].B[obsIndex[t+1]] == 0) {}
				
				probs[t+1][j] = sum * HMM[j].B[obsIndex[t+1]];
			}
		}

		return probs;
	}

	/*
	 * Given an array of hidden Markov states, an array of vocabulary words, an observation set, the number of
	 * words in that set, and the length of an observation sequence, this algorithm calculates the backward
	 * probability of observing the given sequence. It then returns the table of backward probabilities for
	 * each state and time t.   
	 */
	public static double[][] backwardAlgo(State HMM[], String vocab[], String words[], int numWords, int seqLen){

		int numStates = HMM.length;

		// Determine the index of each observation symbol (so we can find it in the A and B matrices)
		int obsIndex[] = new int[seqLen+1];
		for(int i=0; i<numWords; i++){ obsIndex[i+1] = getIndex(words[i], vocab); }

		// probs[t][i] holds the probability of an observation sequence of length <t> starting in state <i>
		double probs[][] = new double[seqLen+1][numStates];

		// Initialize end probabilities to 1 for t = T
		for(int state=0; state<numStates; state++){
			probs[seqLen][state] = 1;
			}

		// Induction Step: t = T-1 ... 1
		for(int t=seqLen-1; t>=0; t--){				// use interval seqLen-2 ... 0 because it is used as an array index
			for(int state=0; state<numStates; state++){

				double sum = 0.0;
				for(int j=0; j<numStates; j++){
					sum += HMM[state].A[j] * HMM[j].B[obsIndex[t+1]] * probs[t+1][j];
				}
				System.out.println("beta: "+ Arrays.deepToString(probs));
				probs[t][state] = sum;
			}
		}

		return probs;
	}

	/*
	 * Given an array of hidden Markov states, an array of vocabulary words, an observation set, the number of
	 * words in that set, and the length of an observation sequence, this algorithm calculates the forward
	 * probability of observing the given sequence. It then returns the table of forward probabilities for
	 * each state and time t.   
	 */
	public static Viterbi viterbiAlgo(State HMM[], String vocab[], String words[], int numWords, int seqLen){

		int numStates = HMM.length;

		// Determine the index of each observation symbol (so we can find it in the A and B matrices)
		int obsIndex[] = new int[numWords+1];
		for(int i=0; i<numWords; i++){ obsIndex[i+1] = getIndex(words[i], vocab); }

		// delta[t][i] holds the probability of seeing observation symbol <t> while in state <i>
		double delta[][] = new double[seqLen+1][numStates];
		int psi[][] = new int[seqLen+1][numStates];

		// Initialize pi probabilities
		for(int state=0; state<numStates; state++){ delta[0][state] = HMM[state].pi; }

		// Initialization Step: t = 1
		for(int state=0; state<numStates; state++){ delta[1][state] = delta[0][state] * HMM[state].B[obsIndex[1]]; }

		// Induction Step: t = 2.. seqLen-1
		for(int t=2; t<seqLen; t++){
			for(int state=0; state<numStates; state++){

				delta[t][state] = findMax(HMM, delta[t-1], state) * HMM[state].B[obsIndex[t]];
				psi[t][state] =  findArgMax(HMM, delta[t-1], state);
			}
		}

		// Termination Step
		double p[] = new double[seqLen+1];		// Holds probability of best path of length t 
		int q[] = new int[seqLen+1];				// Holds state with highest probability at time t

		for(int t=0; t<seqLen; t++){

			p[t] = delta[t][0];
			q[t] = 0;

			for(int i=0; i<numStates; i++){
				if(delta[t][i] > p[t]){
					p[t] = delta[t][i];
					q[t] = i;
				}
			}
		}

		Viterbi v = new Viterbi(p, q);

		return v;
	}

	/*
	 * The getIndex function takes in an observation symbol, and a list of vocabulary words
	 * and return the index of that word in the list. This is used for looking up observation
	 * probabilities of a symbol in a certain State
	 */
	public static int getIndex(String word, String vocab[]){

		int i = 0;
		while(!vocab[i].equalsIgnoreCase(word) && i <= vocab.length) { i++; }

		return i;
	}

	/*
	 * The findMax function takes in an array of States, an array of observation probabilities(for some
	 * specific sequence length), and the index of another state, j. The method determines and returns the 
	 * maximum value of probs[i] * HMM[i].A[j] which is the probability of an observation sequence ending in 
	 * state i multiplied by the probability of transitioning from state i to j (for all i, 1<=i<=N) 
	 */
	public static double findMax(State HMM[], double probs[], int j){

		double max = 0.0;
		for(int i=0; i<HMM.length; i++){

			double x = probs[i] * HMM[i].A[j];
			if( x > max) { max = x; }
		}

		return max;
	}

	/*
	 * The findArgMax function takes in an array of States, an array of observation probabilities(for some
	 * specific sequence length), and the index of another state, j. The method returns the index i
	 * (for all states i, 1<=i<=N) which results in the the maximum value of probs[i] * HMM[i].A[j] which is the
	 * probability of an observation sequence ending in state i multiplied by the probability of transitioning 
	 * from state i to j.
	 */
	public static int findArgMax(State HMM[], double probs[], int j){

		double max = 0;
		int index = 0;
		for(int i=0; i<HMM.length; i++){

			double x = probs[i] * HMM[i].A[j];
			if( x > max) { max = x; index = i; }
		}

		return index;
	}

	/* 
	 * The Viterbi class is used as a wrapper to hold the p and q arrays  
	 */
	public static class Viterbi{

		double p[];		// Probabilities along most likely path
		int q[];		// States that make up most likely path

		public Viterbi(){}
		public Viterbi(double p[], int q[]){
			this.p = p;
			this.q = q;
		}
	}
}