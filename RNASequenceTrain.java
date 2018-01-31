import java.util.*;
import java.io.*;

/** This class implements a Profile Hidden Markov Modell (PHMM) for detection of
	ribosomal RNA-Sequences (rRNA). It is trained with a MSA of trainingseq. and
	calculates the Viterbi-Algortihm for the testsequences.*/
public class RNASequenceTrain {

	
	//enumeration for the states
	private enum State {	
		M(0), I(1), D(2), B(3), E(4);		//Match, Insert, Delete, Begin, End
		public int row;
		State (int row) {this.row = row;}
	}


	//enumeration for the state transitions
	private enum ST {
		MTM(0), MTI(1), MTD(2), ITM(3), ITI(4),		//Match to Match,
			ITD(5), DTM(6), DTI(7), DTD(8);			//Match to Insert, ...
	 	public int row; 									
		ST(int row) { this.row = row; }
		
		public static ST getTransition(State prev, State succ) {
			ST ret = null;
		
			switch (prev) {
				
					case M : //Statetransition M to ..
						switch (succ) {
							
							case M :			//  M
								ret = ST.MTM;
								break;
								
							case I :			//  I
								ret = ST.MTI;
								break;
								
							case D :			//  D
								ret = ST.MTD;
								break;
						
						}
						break;
						
					case I : //Statetransition I to ..
						switch (succ) {
							
							case M :			// M
								ret = ST.ITM;
								break;
								
							case I :			// I
								ret = ST.ITI;
								break;
							
							case D :
								ret = ST.ITD;
								break;
									
						}
						break;
						
					case D : //Statetransition D to ..
						switch (succ) {
							
							case M :			// 	M
								ret = ST.DTM;
								break;
								
							case I :
								ret = ST.DTI;
								break;
								
							case D :			//  D
								ret = ST.DTD;
								break;
							
						}
						break;
						
					default :	//should not occur
						System.err.println("Oh oh, we have a problem!"
							+prev+" "+succ);
						
						if (prev == State.B)	//we abuse MT.. for BT..
							ret = getTransition(State.M, succ);
						else if ( succ == State.E)	//we abused ..TM for ..TE
							ret = getTransition(prev, State.M);
						
						break;
				} 
	
	
			return ret;
		}
		
	}	//end of ST
	
	
	
	//number of statetransitions
	private static final int stsize = ST.values().length;



	//String array for sequence headers
	private static String[] header;
	//String array for sequences
	private static String[] seq;
	//Variable for number of sequences and length of sequence (rows & columns of
	//Matrix)
	private static int seqnr, seqlen;
	
	//boolean array for matchingpositions
	private static boolean[] isMatch;
	//number of matchingpositions
	private static int matnr;
	
	//matrix for counts of emissions at each match- or insertposition
	//ecounts[position][base]
	private static int[][] ecounts;	
	//matrices for estimated emission probabilities (+ pseudocounts) [pos][base]
	//at matching positions, insertion positions & final one for both
	private static double[][] eprobmat, eprobins, eprob;
	
	//matrix for counts of transition at each position [pos][StateTransition]
	private static int[][] tcounts;
	//matrix for estimated transition probabilities (+ pseudoc.) [pos][ST]
	private static double[][] tprob;
	
	//global pseudocount
	private static int pseudo = 1;
		
	
	//matrices for dynamic programming of viterbi variables
	private static double[][] vM, vI, vD;
	//matrices for backtracking
	private static dir[][] btM, btI, btD;
	//constants for directions for backtracking
	private enum dir {
		UPLEFT, UP, LEFT, RIGHT;
	}
	
	//String for the sequence that shall be tested (viterbi)
	private static String testseq;
	//String arrays for testsequences and their names (headers)
	private static String[] tester;
	private static String[] testerheader;
	//number of testsequences
	private static int testernr;
	
	//threshold for logarithmic score
	private static double threshold;
	
	
	
	public static void main(String[] args) throws IllegalArgumentException, 
		IOException {
		
		if (args.length < 1) {
			System.out.println("The command line parameters must be at least "
				+"the path to the file with the MSA for the training of the HMM"
				+" and optional a path to a file with testsequences!");
			System.exit(1);
		}
		
		
		
		String path = args[0];
		String line = null;
		ArrayList<String> array = new ArrayList<String>();
		BufferedReader r = new BufferedReader(new FileReader(path));
		
		//read data from the file
		while ( (line = r.readLine()) != null )
			if (!line.startsWith(";")) //ignore comment lines
				array.add(line);
		
		r.close();
		
		
		//init arrays for headers and sequences
		seqnr = array.size() / 2;
		header = new String[seqnr];
		seq = new String[seqnr];
		
		//split up the input into the right array
		for (int i = 0; i < header.length; i++)
			header[i] = array.get(2*i);
		for (int i = 0; i < seq.length; i++)
			seq[i] = array.get(2*i+1);
			
		
		//train the PHMM & print the emissionprobabilities for matchpositions
		train();
		
		if (args.length == 1) {
			printMatchPositions();
			System.exit(0);		//end here
		}	
		
		
		
		
		//else test sequences from path in args[1]
		System.out.println("Starting to test...");
		path = args[1];
		line = null;
		array = new ArrayList<String>();
		r = new BufferedReader(new FileReader(path));
		
		while ( (line = r.readLine()) != null)
			if (!line.startsWith(";"))
				array.add(line);
		
		r.close();
		
		//init tester & testerheader
		testernr = array.size() / 2;
		tester = new String[testernr];
		testerheader = new String[testernr];
		
		for (int i = 0; i < testerheader.length; i++)
			testerheader[i] = array.get(2*i);
		for (int i = 0; i < tester.length; i++)
			tester[i] = array.get(2*i+1);
		
		
		//now construct a threshold for accepting/denying the testsequences with
		//the trainingsequences
		//threshold = minimumThreshold(5, 0.9);
		
		
		//test the testsequences and decide rRNA/Non-rRNA
		System.out.println("Testingsequences:");
		double temp;
		for (String t : tester) {
			testseq = t;
			temp = processViterbi();			
				
			//for matchhitThreshold
			boolean isrRNA = matchhitThreshold(0.8, 3);
			System.out.println(temp+"\t"+ (isrRNA ? 1 : 0));
		}
	
	
	}
		
		
	
	/** Prints the matchpositions and the emissionprobabilities for each char.*/
	private static void printMatchPositions() {
		
		System.out.println("The matchpositions with the character and the "
				+ "emissionprobabilities for this matchposition");
		int ptr = 0;
		System.out.println("Pos.\te_pos('A')\te_pos('C')\te_pos('G')"
						+"\te_pos('U')");
		for (int i = 0; i < isMatch.length; i++) {
			if (isMatch[i]) {
				System.out.print((ptr+1)+"\t");
				
				for (int l = 0; l < eprobmat[0].length; l++)
					System.out.format("%.8f\t", eprobmat[ptr][l]);
					
				System.out.println();
				ptr++;
			}
		}
		
		System.out.println();
	
	}
		
		
		
		
	
	/** This method trains the profile hidden markov model:
		1. Get the structure of the PHMM.
		2. Count the emissions in the MSA & estimate the emissionprobabilities.
		3. Count the transitions in the MSA & estimate the transitionprobs.
		*/	
	private static void train() {
		
		System.out.println("#Starting to train......");
		
		
		//1. Get matchpositions (more than 50% gaps in one column -> no matchp.)
		System.out.println("Getting matchpositions...");
		seqlen = seq[0].length();
		isMatch = new boolean[seqlen];
		matnr = 0;
		int lastmatchpos = 0;
		for (int j = 0; j < seqlen; j++) {
			int countgaps = 0;
			for (int i = 0; i < seqnr; i++){
				if (seq[i].charAt(j) == '-')
					countgaps++; 	
			}
			if (2 * countgaps >= seqnr)	//matchposition?
				isMatch[j] = false;			//no
			else {
				isMatch[j] = true;			//yes
				matnr++;
				lastmatchpos = j;
			}
			
		}
		
		

		System.out.println("Got matchpositions ("+matnr+")!");

		//2. Count all emissions for each pos. & estimate emissionprobabilities
		System.out.println("Counting all emissions...");
		ecounts = new int[seqlen][4];
		eprobmat = new double[matnr][4];
		eprobins = new double[matnr+1][4]; 	//before and after every matchpos.
		
		
		//2.1 Count emissions for each column
		for (int j = 0; j < seqlen; j++) {
			int[] temp = new int[4];   //[0]='A',[1]='C',[2]='G',[3]='U' 
			for (int i = 0; i < seqnr; i++)
				switch (seq[i].charAt(j)) {
					case 'A' : 
						temp[0]++;
						break;
						
					case 'C' :
						temp[1]++;
						break;
						
					case 'G' :
						temp[2]++;
						break;
						
					case 'U' : 
						temp[3]++;
						break;
						
					case '-' :
						//temp[4]++;	//ignore gaps		
						break;
						
					default :
						throw new IllegalArgumentException("The character "+
						seq[i].charAt(j)+" in sequence "+i+" at position "+j+
						" is not valid");
				}
				
			ecounts[j] = temp;
				 
		}
		
		
		System.out.println("Counted all emissions!");
		
		
		//2.2 Estimate emissionprobabilities with relative frequencies
		//2.2.1 Estimate at matchingpositions
		System.out.println("Estimating at matchpositions...");
		int index = 0;
		for (int j = 0; j < eprobmat.length; j++) {
			int sum = 0;
			
			//get the next matching position
			while (!isMatch[index])
				index++;
			
			//Matchingposition
			//get the sum of all counts on this column
			for (int i = 0; i < ecounts[index].length; i++)
				sum += ecounts[index][i];
		
			//estimate probabilities with laplace-rule
			double[] p = new double[ecounts[index].length];
			for (int i = 0; i < p.length; i++)
				p[i] = (double) (ecounts[index][i] + pseudo)
						/ (sum + p.length * pseudo);
			
			//set the column in the matrix & increment index
			eprobmat[j] = p;
			index++;
		}
		
		
		System.out.println("Estimated at matchpositions!");
		

		//2.2.2 Estimate at insertionpositions
		System.out.println("Estimating at insertpositions...");
		index = 0;
		
		boolean flag = false;
		for (int j = 0; j < eprobins.length; j++) {
			int sum = 0;
			int end = 0;

			//find next matchingposition
			end = index + 1;
			if (flag || index > 0)
				while (end < isMatch.length && !isMatch[end])
					end++;
			else
				if (isMatch[0]) {
					index = -1;		
					end = 0;
					flag = true;
				}
				else
					while (end < isMatch.length && !isMatch[end])
						end++;
			
			
			
			if (end > index + 1) {
				//get sum of all ins.positions (columns) from index to end - 1
				for (int pos = index; pos < end; pos++) {
					for (int i = 0; i < ecounts[pos].length; i++)
						sum += ecounts[pos][i];	
				}
				
				//make one column by adding all together
				long[] col = new long[ecounts[index].length];
				for (int pos = index; pos < end; pos++) {
					for (int i = 0; i < ecounts[pos].length; i++)
						col[i] += ecounts[pos][i];		
				}
				
				//estimate emissionprobabilities with laplace-rule
				double[] p = new double[ecounts[index].length];
				for (int i = 0; i < p.length; i++)
					p[i] = (double) (col[i] + pseudo) 
							/ (sum + p.length * pseudo);
					
				eprobins[j] = p;
			}
			else {
				//one insertpos. betw. two matchpos. -> no emissions only pseudo
				double[] p = new double[ecounts[0].length];
				for (int i = 0; i < p.length; i++)
					p[i] = (0.0 + pseudo) 
							/ (sum + p.length * pseudo);
				eprobins[j] = p;
			}
			
			
			index = end;
		}

		
		System.out.println("Estimated at insertpositions!");
		
		
		//save emissionprobabilities in one matrix, first matchemissions,
		//second insertemissions
		eprob = 
			new double[eprobins.length][eprobmat[0].length+eprobins[0].length];
		
		for (int i = 0; i < eprob.length; i++) {
			
			//first matchemissions
			if (i > 0)
				for (int j = 0; j < eprobmat[0].length; j++) 
					eprob[i][j] = eprobmat[i-1][j];
			
			//second insertemissions
			for (int j = 0; j < eprobins[0].length; j++)
				eprob[i][j+eprobmat[0].length] = eprobins[i][j];
		
		}
		
		System.out.println("Saved emissionprobabilities to one matrix!");
		
		
		
		
		//2.3 Count all transitions & estimate transitionprobabilities
		System.out.println("Counting all transitions...");
		tcounts = new int[matnr + 1][stsize];	//plus one for initial insert
		
		for (int i = 0; i < seqnr; i++) {
			
			//2.3.1 Build state path for this sequence
			State[] p = new State[seqlen];
						
			for (int j = 0; j < p.length; j++) {
				if (isMatch[j]) {	//it is a matchposition
					if (seq[i].charAt(j) == '-')//gap means deletion
						p[j] = State.D;
					else						//no gap means matching
						p[j] = State.M;			//must match because MSA
				}
				else {				//it is no matchposition
					if (seq[i].charAt(j) == '-')
						p[j] = null;			//gap will be skipped
					else
						p[j] = State.I;			//no gap means insertion
				}
			}
		
		
			//2.3.2 Count transitions for this path
			
			
			State succ = null, s = p[0];
			int ptr = 0;
			
			//first column of tcounts are transitions from State.B to[0]			
			if (s == State.M)	 			//abuse MTM row for BTM
				tcounts[0][ST.MTM.row]++;
			else if (s == State.I)			//abuse MTI row for BTI 
				tcounts[0][ST.MTI.row]++;		//(!ITI is needed)
			else if (s == State.D)			//abuse MTD row for BTD
				tcounts[0][ST.MTD.row]++;	 			
			else if (s == null) {
				//were starting with null so get the first not null state
				int l = 1;		//we shouldnt have a seq like "----.....-----"
				while (p[l] == null)
					l++;
				s = p[l];
				if (s == State.M)		//look above
					tcounts[0][ST.MTM.row]++;
				else if (s == State.I)			
					tcounts[0][ST.MTI.row]++;		
				else if (s == State.D)		
					tcounts[0][ST.MTD.row]++;
				else {
					System.out.println("Problem at tcounts[0]!");
					
				}
			}
			
			//is first state an insertstate with gap(s)?
			if (p[0] == null) {
				//find the succ. not null state and check it
				int l = 1;
				while (p[l] == null)
					l++;
				s = p[l];
				
				if (s == State.M || s == State.D)
					ptr = 1;
				else
					ptr = 0;
			}
			
			for (int j = 0; j < p.length - 1; j++) {
				s = p[j];
				
				//skip because there was no transition
				if (s == null)
					continue;
				
				
				//is first state a match or delete?
				if (j == 0) {
					if (s == State.M || s == State.D)
						ptr = 1;
				}
				
				//get succ. state
				int k = j+1;
				while (k < p.length && p[k] == null)
					k++;
				
				if (k >= p.length) {
					//path ends with a (chain of) null value(s)
					//will be handled after for loop
					break;
				}	
		 	
		 		succ = p[k];
					 
				ST st = ST.getTransition(s, succ);
				if (st == null) {	
					//should not occur
					System.err.println("NULL at ST.getTransition(...)!");
					System.exit(1);
				}
				
				
				tcounts[ptr][st.row]++;
				
				
				
			
				//only increment ptr if the transition was not to an insertstate
				if (st == ST.MTI || st == ST.ITI || st == ST.DTI)
					;
				else
					ptr++;
			}
			
			
			//last column are transitions to end state
			//tcounts[seqlen][ST.getTransition(p[p.length-1], State.E).row]++;
			s = p[p.length - 1];
			
			if (s == State.M)			//here abuse: the ..TM states
				tcounts[tcounts.length - 1][ST.MTM.row]++;
			else if (s == State.D)
				tcounts[tcounts.length - 1][ST.DTM.row]++;
			else if (s == State.I)
				tcounts[tcounts.length - 1][ST.ITM.row]++;
			else if (s == null) {
				//this is the case with null(chain) at the end; analog to start:
				//find the first prev. not null state -> count transition to end
				int l = p.length - 2;
				while (p[l] == null)
					l--;
				s = p[l];
				
				if (s == State.M)			//see above
					tcounts[tcounts.length - 1][ST.MTM.row]++;
				else if (s == State.D)
					tcounts[tcounts.length - 1][ST.DTM.row]++;
				else if (s == State.I)
					tcounts[tcounts.length - 1][ST.ITM.row]++;
			
			}
		
		}	//end of for over all paths
		
		
		
		
		
		
		//for (ST s : ST.values())
			
		
		
		
		System.out.println("Counted all transitions!");
		
		
		
		
		//2.3.3 Estimate the transitionprobabilities out of the counts
		System.out.println("Estimating transitionprobabilities...");
		tprob = new double[tcounts.length][];
		
		for (int i = 0; i < tcounts.length - 1; i++) {
			double[] p = new double[tcounts[0].length];
			int sum = 0; int j;
			
			//first sum over all from match outgoing transitions -> estimate
			for (j = ST.MTM.row; j <= ST.MTD.row; j++) 
				sum += tcounts[i][j];
			
			for (j = ST.MTM.row; j <= ST.MTD.row; j++)
				p[j] = ((double) tcounts[i][j] + pseudo) 
							/ (sum + (ST.MTD.row - ST.MTM.row + 1) * pseudo); 
				
				
			//second sum over all from insert outgoing transitions
			sum = 0;
			for (j = ST.ITM.row; j <= ST.ITD.row; j++) 
				sum += tcounts[i][j];
				
			for (j = ST.ITM.row; j <= ST.ITD.row; j++)
				p[j] = ((double) tcounts[i][j] + pseudo) 
							/ (sum + (ST.ITD.row - ST.ITM.row + 1) * pseudo);
				
				
			//third sum over all from delete outgoing transitions
			sum = 0;
			for (j = ST.DTM.row; j <= ST.DTD.row; j++)
				sum += tcounts[i][j];
			
			for (j = ST.DTM.row; j <= ST.DTD.row; j++)
				p[j] = ((double) tcounts[i][j] + pseudo) 
							/ (sum + (ST.DTD.row - ST.DTM.row + 1) * pseudo);
			 
		
			tprob[i] = p;
		
		}
	
	
		//last column needs special treatment -> ..TM transitions are ..TE here
		int lastcol = tcounts.length - 1;
		tprob[lastcol] = new double[tprob[0].length];
		int sum = tcounts[lastcol][ST.MTM.row];
		sum += tcounts[lastcol][ST.MTI.row];
		
		tprob[lastcol][ST.MTM.row] = ((double) tcounts[lastcol][ST.MTM.row] 
			+ pseudo) / (sum + 2 * pseudo);
		tprob[lastcol][ST.MTI.row] = ((double) tcounts[lastcol][ST.MTI.row] 
			+ pseudo) / (sum + 2 * pseudo);
		
		sum = tcounts[lastcol][ST.ITM.row];
		sum += tcounts[lastcol][ST.ITI.row];
		
		tprob[lastcol][ST.ITM.row] = ((double) tcounts[lastcol][ST.ITM.row] 
			+ pseudo) / (sum + 2 * pseudo);
		tprob[lastcol][ST.ITI.row] = ((double) tcounts[lastcol][ST.ITI.row] 
			+ pseudo) / (sum + 2 * pseudo);
		
		sum = tcounts[lastcol][ST.DTM.row];
		sum += tcounts[lastcol][ST.DTI.row];
		
		tprob[lastcol][ST.DTM.row] = ((double) tcounts[lastcol][ST.DTM.row] 
			+ pseudo) / (sum + 2 * pseudo);
		tprob[lastcol][ST.DTI.row] = ((double) tcounts[lastcol][ST.DTI.row] 
			+ pseudo) / (sum + 2 * pseudo);
			
		
		System.out.println("Estimated transitionprobabilities!");
		System.out.println("#Finished training!\n");
		
	
	}
	

	/** Calculates vM at the position (i, j).*/
	private static void viterbiM(int i, int j) {
		double maximum = - Double.MAX_VALUE;
		double[] potMax = new double[3];
		int maxCase = 0;
		
		
		//calculate the cases
		potMax[0] = vM[i-1][j-1] + Math.log(tprob[j-1][ST.MTM.row]);
		potMax[1] = vI[i-1][j-1] + Math.log(tprob[j-1][ST.ITM.row]);
		potMax[2] = vD[i-1][j-1] + Math.log(tprob[j-1][ST.DTM.row]);
		
	
	
		//get the maximum
		for (int z = 0; z < potMax.length; z++)
			if (maximum < potMax[z]) {
				maximum = potMax[z];
				maxCase = z;
			}
		
		//multiplicate with emissionprobability
		int xi = "ACGU".indexOf(testseq.charAt(i - 1));
		vM[i][j] = Math.log(eprobmat[j - 1][xi]) + maximum;
		
		//for debug/backtrack
		dir d = null;
		switch (maxCase) {
			case 0 : d = dir.UPLEFT; break;
			case 1 : d = dir.UP; break;
			case 2 : d = dir.LEFT; break;
		}
		btM[i][j] = d;
	}
	
	
	
	/** Calculates vI at the position (i, j).*/
	private static void viterbiI(int i, int j) {
		double maximum = - Double.MAX_VALUE;
		double[] potMax = new double[3];
		int maxCase = 0;
		
		//calculate the cases
		potMax[0] = vM[i-1][j] + Math.log(tprob[j][ST.MTI.row]);
		potMax[1] = vI[i-1][j] + Math.log(tprob[j][ST.ITI.row]);
		potMax[2] = vD[i-1][j] + Math.log(tprob[j][ST.DTI.row]);
		
	
	
		//get the maximum
		for (int z = 0; z < potMax.length; z++)
			if (maximum < potMax[z]) {
				maximum = potMax[z];
				maxCase = z;
			}
		
		//multiplicate with emissionprobability
		int xi = "ACGU".indexOf(testseq.charAt(i - 1));
		vI[i][j] = Math.log(eprobins[j][xi]) + maximum;
		
		//for debug/backtrack 
		dir d = null;
		switch (maxCase) {
			case 0 : d = dir.UPLEFT; break;
			case 1 : d = dir.UP; break;
			case 2 : d = dir.LEFT; break;
		}
		btI[i][j] = d;
	}
	
	
	
	/** Calculates vD at the position (i, j).*/
	private static void viterbiD(int i, int j) {
		double maximum = - Double.MAX_VALUE;
		double[] potMax = new double[3];
		int maxCase = 0;
	
		//calculate the cases
		potMax[0] = vM[i][j-1] + Math.log(tprob[j-1][ST.MTD.row]);
		potMax[1] = vI[i][j-1] + Math.log(tprob[j-1][ST.ITD.row]);
		potMax[2] = vD[i][j-1] + Math.log(tprob[j-1][ST.DTD.row]);
		
		
		//get the maximum
		for (int z = 0; z < potMax.length; z++)
			if (maximum < potMax[z]) {
				maximum = potMax[z];
				maxCase = z;
			}
		
		//delete is 'silent' -> no emissionprobability
		vD[i][j] = maximum;
		
		//for debug/backtrack
		dir d = null;
		switch (maxCase) {
			case 0 : d = dir.UPLEFT; break;
			case 1 : d = dir.UP; break;
			case 2 : d = dir.RIGHT; break;
		}
		btD[i][j] = d;
	}
	

	/** Process the viterbi algorithm for the testsequence in the variable
		testseq.
		@return Score of the viterbi path*/
	private static double processViterbi() {
		//init v's & bt's and calculate the viterbipath with them
		int lengthOfSequence = testseq.length();
		vM = new double[1 + lengthOfSequence][1 + matnr];
		btM = new dir[vM.length][vM[0].length];
		vI = new double[1 + lengthOfSequence][1 + matnr];
		btI = new dir[vI.length][vI[0].length];
		vD = new double[1 + lengthOfSequence][1 + matnr];
		btD = new dir[vD.length][vD[0].length];
		
		
		//starting here log(1.0) = 0, first row and col -inf
		vM[0][0] = 0; 	
		for (int i = 1; i < vM.length; i++)
			vM[i][0] = Double.NEGATIVE_INFINITY;
		for (int j = 1; j < vM[0].length; j++)
			vM[0][j] = Double.NEGATIVE_INFINITY;
			
		for (int j = 0; j < vI.length; j++) 
			vI[j][0] = Double.NEGATIVE_INFINITY;
		
		for (int i = 0; i < vD[0].length; i++)
			vD[0][i] = Double.NEGATIVE_INFINITY;
		
		
		
		for (int i = 1; i < vM.length; i++) {
			//viterbiD(0, 1);
			for (int j = 1; j < vM[0].length; j++) {
				viterbiM(i, j);
				viterbiI(i, j);
				viterbiD(i, j);
			}
		}
		
		
		
	
	
		//termination
		State ter;
		double max;
		if (vM[lengthOfSequence][matnr] >= vI[lengthOfSequence][matnr] &&
			vM[lengthOfSequence][matnr] >= vD[lengthOfSequence][matnr]) {
			ter = State.M;
			max = vM[lengthOfSequence][matnr];
		
		} else if (vI[lengthOfSequence][matnr] >= vD[lengthOfSequence][matnr]) {
			ter = State.I;
			max = vM[lengthOfSequence][matnr];
		} else {
			ter = State.D; 	
			max = vM[lengthOfSequence][matnr];
		}
		
		//Debug
		
		
		return max;
		}
	
	/** Calculate the minimum score of the first n trainingsequences plus the
		logarithmized factor
		@param n The first n trainingsequence's scores will be compared
		@param factor Factor for threshold
		@return logar. Threshold (incl. factor)*/
	private static double minimumThreshold(int n, double factor) {
		
		double min = Double.MAX_VALUE, temp;
	
		for (int i = 0; i < n; i++) {
			testseq = seq[i].replace("-","");
			temp = processViterbi();
			if (min > temp)
				min = temp;
		}
		
		return min + Math.log(factor);
	}
	
	/** Decide, whether the testseq is rRNA / Non-rRNA (see protocol) with the
		"matchhit-method".
		@param qthr Threshold for Q
		@param lthr Threshold for L
		@return True, if (Q >= qthr && L >= lthr) --> rRNA*/
	private static boolean matchhitThreshold(double qthr, double lthr) {
			
		ArrayList<State> path = new ArrayList<State>();
		State ter;
		int lengthOfSequence = testseq.length();
		//get terminating state
		if (vM[lengthOfSequence][matnr] >= vI[lengthOfSequence][matnr] &&
			vM[lengthOfSequence][matnr] >= vD[lengthOfSequence][matnr])
			ter = State.M;
		else if (vI[lengthOfSequence][matnr] >= vD[lengthOfSequence][matnr])
			ter = State.I;
		else 
			ter = State.D; 	
		
		path = backtrack(ter, lengthOfSequence, matnr, path);
		
	
		//count matchhits of path
		int hitcounter = 0;
		for (int i = 0; i < path.size(); i++)
			if (path.get(i) == State.M)
				hitcounter++;
	
		//count chainlengths
		ArrayList<Integer> chainlengths = new ArrayList<Integer>();
		int index = 0, end;
		while (index < path.size()) {
			if (path.get(index) != State.M) { 	//no match so go on
				index++;	
			}
			else {
				//find end of chain
				end = index + 1;
				while (end < path.size() && path.get(end) == State.M)
					end++;
				
				chainlengths.add(end - index);
				
				index = end + 1;
			}
		
		}
		
		//calculate Q and L
		double Q = (double) hitcounter / lengthOfSequence;
		double L = 0;
		for (int i = 0; i < chainlengths.size(); i++)
			L += chainlengths.get(i);
		L /= chainlengths.size();
		
		
		
		return (Q >= qthr && L >= lthr);
	}
	
	
	
	
	/** Backtrack the viterbipath recursively.
		@param ter Terminal state
		@param i Character at position i in testseq
		@param j State at position j
		@param list ArrayList where the path shall be saved in
		@return Viterbipath*/
	private static ArrayList<State> backtrack(State ter, int i, int j,
				ArrayList<State> list) {
		if (i == 0 && j == 0)
			return list;
	
		if (i == 0) {
			//there are no more characters -> so there are deletions
			while (j > 0) {
				list.add(State.D);
				j--;
			}
			return list;
		}
		else if (j == 0) {
			//there are still characters -> so there are insertions
			while (i > 0) {
				list.add(State.I);
				i--;
			}
			return list;
		}
	
	
	
		State st = null;
		//get the direction we came from -> recursion
		switch (ter) {
			
			case M :
				if (btM[i][j] == dir.UPLEFT)
					st = State.M;
				else if (btM[i][j] == dir.UP)
					st = State.I;
				else
					st = State.D;
					
				list = backtrack(st, i-1, j-1, list);
				break;
				
			case I :
				if (btI[i][j] == dir.UPLEFT)
					st = State.M;
				else if (btI[i][j] == dir.UP)
					st = State.I;
				else 
					st = State.D;
					
				list = backtrack(st, i-1, j, list);
				break;
				
			case D : 
				if (btD[i][j] == dir.UPLEFT)
					st = State.M;
				else if (btD[i][j] == dir.RIGHT)
					st = State.D;
				else 
					st = State.I;
					
				list = backtrack(st, i, j-1, list);
				break;
			
		}
		
		//return path plus the last state (ter)
		
		list.add(ter);
		return list;
	
	}
	
	
	
	
	
	
	
	/** For debug/backtrack purposes. Returns the state path to the last viterbi
		process with backtracking recursively.
		@return State path*/
	private static String printStateSeq(State ter, int i, int j) {
		if (i == 0 && j == 0)
			return "";
		
		if (i == 0) {
			StringBuilder sb = new StringBuilder();
			while (j > 0) {
				sb.append(j+"D ");
				j--;	
			}
			return sb.reverse().toString();
		}
			
		if (j == 0) {
			StringBuilder sb = new StringBuilder();
			while (i > 0) {
				sb.append("I0 ");
				i--;
			}
			return sb.toString();
		}
			
		
		String ret = "";
		State st = null;
	
		switch (ter) {
			
			case M :
				if (btM[i][j] == dir.UPLEFT)
					st = State.M;
				else if (btM[i][j] == dir.UP)
					st = State.I;
				else
					st = State.D;
					
				ret = printStateSeq(st, i-1, j-1) + " " + ter + j +ret;
				break;
				
			case I :
				if (btI[i][j] == dir.UPLEFT)
					st = State.M;
				else if (btI[i][j] == dir.UP)
					st = State.I;
				else 
					st = State.D;
					
				ret = printStateSeq(st, i-1, j) + " " + ter + j + ret;
				break;
				
			case D : 
				if (btD[i][j] == dir.UPLEFT)
					st = State.M;
				else if (btD[i][j] == dir.RIGHT)
					st = State.D;
				else 
					st = State.I;
					
				ret = printStateSeq(st, i, j-1) + " " + ter + j + ret;
				break;
			
		}
	
		return ret;
	
	} 
	
	
	
}
