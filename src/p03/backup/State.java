package p03.backup;

/*
 * State.java: Defines a State in the Hidden Markov Model
 */

public class State {

	String name;		// The name of the state 
	int num;			// State number
	double A[];			// The A matrix hold the transition probabilities to each state from this one
	double B[];			// The B matrix hold the probabilities of observing each vocab word in this state
	double pi;
	
	public State(){
		
	}
	
	public State(String n, int x, String a, String b, double p){
		name = n;
		num = x;

		// Split a into tokens
		String splitA[] = a.split(" ");
		A = new double[splitA.length];
		
		// Fill A array for this state
		for(int i=0; i<splitA.length; i++){
			A[i] = Double.parseDouble(splitA[i]);
		}
		
		// Split b into tokens
		String splitB[] = b.split(" ");
		B = new double[splitB.length];
		
		// Fill B array for this state
		for(int j=0; j<splitB.length; j++){
			B[j] = Double.parseDouble(splitB[j]);
		}
		
		pi = p;		// Set probability of starting in this state 
	}
	
	// Alter the transition probability from this state
	public void changeA(int state, double value){
		A[state] = value;
	}
	
	// Alter the probability of seeing a vocabulary word in this state
	public void changeB(int word, double value){
		A[word] = value;
	}
	
	// Alter the probability of starting in this state
	public void changePi(double value){
		pi = value;
	}
	
}
