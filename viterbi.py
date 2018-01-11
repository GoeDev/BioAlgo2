import numpy as np

'''
N: number of hidden states
'''
class Decoder(object):
    #Konstruktor für Decoder-Klasse (mittels viterbi)
    def __init__(self, startWahr, transWahr, beobWahr):
        self.N = startWahr.shape[0]
        self.startWahr = startWahr
        self.transWahr = transWahr
        self.beobWahr = beobWahr
        assert self.startWahr.shape == (self.N, 1)
        assert self.transWahr.shape == (self.N, self.N)
        assert self.beobWahr.shape[0] == self.N
        
    #gibt die beobachtete Wahrscheinlichkeit in der gewählten Spalte zurück
    def Beob(self, beob):
        return self.beobWahr[:, beob, None]

    def Decode(self, beob):
        gitter = np.zeros((self.N, len(beob)))
        backpt = np.ones((self.N, len(beob)), 'int32') * -1
                
        #Gitter initialisieren
        gitter[:, 0] = np.squeeze(self.startWahr * self.Beob(beob[0]))
                
        for t in xrange(1, len(beob)):
            gitter[:, t] = (gitter[:, t-1, None].dot(self.Beob(beob[t]).T) * self.transWahr).max(0)
            backpt[:, t] = (np.tile(gitter[:, t-1, None], [1, self.N]) * self.transWahr).argmax(0)
        #Terminierung
        tokens = [gitter[:, -1].argmax()]
        for i in xrange(len(beob)-1, 0, -1):
            tokens.append(backpt[tokens[-1], i])
        return tokens[::-1]
