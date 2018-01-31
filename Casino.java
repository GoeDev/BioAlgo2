import java.util.*;
import java.io.*;


/** This class implements the Viterbi-Algorithm for the HMM of the 
	'occasionally dishonest casino' (Example Durbin et al. "Biological sequence 
	analysis").*/
public class Casino {

	//Hidden states (F = fair, L = loaded)
	private static char[] hidden_states = {'F','L'};
	
	//State transition matrix with start probabilities, [state1][state2]
	private static double[][] transition_prob = 
			{{0, 0.5,  0.5 },  	//start to start, fair, loaded
			 {0, 0.95, 0.05},	//fair to    -     "      -
			 {0, 0.1 , 0.9 }};	//loaded to  -     "      -
		
	//Emission probability matrix 
	private static double[][] emission_prob = 
		{{1.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6},	//fair 
		 {0.1,   0.1,   0.1,   0.1,   0.1,   0.5  }} ;  //loaded
		 
			
	//Observations, Char [] -> String
	private static String obs;

	//Matrix with all viterbi variables in it (v[i][j] -> viterbi variable at 
	//position j and state i
	private static double[][] v;
	
	//Matrix for backtracking
	private static char[][] bt;

	
	/** Calculate the viterbi variable of a given state at a given position in
		the observation sequence.
		@param statenr The number of the state (0 = Fair, 1 = Loaded (Array))
		@param pos The position in the observation sequence (i-th character)*/
	private static void viterbi(int statenr, int pos) {
		double maximum = - Double.MAX_VALUE, temp;
		char prev_state = 'X';
		
		//for all states get product of viterbi variable up to pos-1 and the
		//transition from that state to the next state (let out start)
		if (pos > 1) {
			for (int i = 1; i < v.length; i++) {
				temp = v[i][pos-1] + Math.log(transition_prob[i][statenr + 1]);  
				if (maximum < temp) {
					maximum = temp;
					prev_state = hidden_states[i-1]; //i - 1 = no start state
				}
			
			}
		}	
		else { //it is the second column (first one with observations) (start)
			maximum = /*log (1.0)+*/ Math.log( transition_prob[0][statenr + 1]);
		}
		
		
		v[statenr + 1][pos] = 
			Math.log(emission_prob[statenr][obs.charAt(pos-1) - '1']) + maximum;
		bt[statenr + 1][pos] = prev_state;
	}

	
	/** Traceback the viterbi path with help of the backtrack matrix bt.
		@return A string containing the traceback of the viterbi path.*/
	private static String traceback() {
		String res = ""; 
		int maxpos = 1;
		int actstate = 0;
		
		//get maxpos
		for (int i = 2; i < v.length; i++)
			if (v[maxpos][v[0].length - 1] < v[i][v[0].length - 1])
				maxpos = i;
			
		res += hidden_states[maxpos-1]; 	//add last state
		actstate = maxpos;
		for (int i = bt[0].length - 1; i > 1; i--) {
			char c = bt[actstate][i];
			res += c;
			actstate = (new String(hidden_states)).indexOf(c)+1;
		}
		
		return (new StringBuilder(res)).reverse().toString();
	}


	/** Main method.*/
	public static void main(String[] args) throws IOException{
		
		if (args.length != 1) 
			throw new IllegalArgumentException("The command line parameter must"
				+" be the path to file with the observations in it!.");
		
		//args[0] should contain path to file with observations
		String path = args[0];
		
		//load observations from file
		String line;
		StringBuilder b = new StringBuilder();
		BufferedReader r = new BufferedReader(new FileReader(path));
		while ((line = r.readLine()) != null)
			b.append(line);
		r.close();
		obs = b.toString();	
		
		v = new double[hidden_states.length + 1][obs.length() + 1];
		bt = new char[hidden_states.length + 1][obs.length() + 1];
		
		//init viterbi matrix
		v[0][0] = 1; //and v[k][0] = 0 (default) for k > 0
		
		//compute all the viterbi variables iteratively
		for (int i = 1; i < v[0].length; i++)
			for (int j = 0; j < hidden_states.length; j++)
				viterbi(j ,i);

		
		
		
		//traceback
		String pi = traceback();
		int l = 60; int start = 0;
		//print the result in lines of length 60
		while (obs.length() - (start + l) >= 0) {
			System.out.println(obs.substring(start, start+l)+"\n"
				+pi.substring(start, start+l)+"\n");
			start += l;
		}
		
		//if there are some characters missing print them
		if (obs.length() > start)
			System.out.println(obs.substring(start)+"\n"+pi.substring(start));
	}	

}
