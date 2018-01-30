import argparse
import math
import sys

#hidden states, F = fair, L = loaded
hiddenStates = ['F', 'L']
#state transition matrix (probabilities)
stateTransitionProb = [[0,0.5,0.5],[0,0.95,0.05],[0,0.1,0.9]]
#emission probability matrix
emissionProb = [[1.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6],[0.1, 0.1, 0.1, 0.1, 0.1, 0.5]]
#string for observations
obs = ""

#calculate viterbi variable of given state and position in the sequence
#statenr = state mapped to number (0=F, 1=L)
#pos = position in the observation sequence
def viterbi():
        maximum = -3.0
        v = [[0. for i in range(len(obs) + 1)] for j in range(len(hiddenStates) + 1)]
        bt = [['' for i in range(len(obs) + 1)] for j in range(len(hiddenStates) + 1)]
        v[0][0]=1
        #iteratively calculate all viterbi variables
        for pos in range(1, len(obs)+1):
                for statenr in range(0,len(hiddenStates)):
                        previousState = 'X' #just init with whatever

                        #iterate over states, get product of viterbi variable until pos-1
                        #and state transitions

                        if pos > 1:
                                for i in range(1, len(v)):
                                        temp = v[i][pos-1] + math.log(stateTransitionProb[i][statenr + 1])
                                        print("temp: " + str(temp) + ", max: " + str(maximum))
                                        if maximum < temp:
                                                maximum = temp
                                                previousState = hiddenStates[i-1]
                        #else: #it has to be the second column
                        #       maximum = math.log(stateTransitionProb[0][statenr + 1])

                        #fill matrices
                        index = statenr + 1
                        ####DEBUG
                        #print("state:"+str(statenr)+", pos:"+str(pos))
                        v[index][pos] = math.log(emissionProb[statenr][int(obs[pos-1])-1]) + maximum
                        bt[index][pos] = previousState

        return v, bt


#method for tracing back the viterbi path with traceback matrix back
#returns string containing the traceback
def traceback(vit, back):
        #initialize variables
        ret = ""
        maxpos = 1
        currentState = 0
        #find out max position
        for i in range(2, len(vit)):
                if vit[maxpos][len(obs) - 1] < vit[i][len(obs) - 1]:
                        maxpos = i

        #add last state
        ret = ret + hiddenStates[maxpos - 1]

        currentState = maxpos
        #count down to 2 from backtracking length
        for i in reversed(range(2, len(back[0]))):
                c = back[currentState][i]
                ret = ret + c
                currentState = hiddenStates.index(c) + 1

        #return reversed string
        return ret[::-1]


#####Main funktion starts here#####



#parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("filename", help="Path to the file containing the sequences")
args = parser.parse_args()



#read input file
with open(args.filename, 'r') as myfile:
        obs=myfile.read().replace('\r\n', '').replace('\n', '')

v, bt = viterbi()
#calculate traceback
pi = traceback(v, bt)

l=60
start=0

while len(obs) - (start + 1) >=0:
        print(obs[start:start+1]+ "\n" + pi[start:start+1] + "\n")
        start = start + 1

#if neccessary, print missing characters
if len(obs) > start:
        print(obs[start] + "\n" + pi[start])
