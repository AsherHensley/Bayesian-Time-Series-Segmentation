# Nonparametric-Learning-for-HMMs-with-Preferential-Attachment

This work considers the idea of Hidden Markov Models where the state transition probabilities change according to a preferential attachment process called the Yule-Simon process. It's assumed there are an infinite number of hidden states, and each time there is a state transition a new state is randomly created. The model is nonparametric in this sense because the number of states is theoretically infinite and is determined by the data during inference. The code in this repositiory creates a synthetic data set and applies a Gibbs sampling procedure to infer the hidden state transition points and Yule-Simon parameter as seen in the paper:

Hensley, Asher A., and Petar M. Djurić. "Nonparametric learning for Hidden Markov Models with preferential attachment dynamics." Acoustics, Speech and Signal Processing (ICASSP), 2017 IEEE International Conference on. IEEE, 2017.

To get started, run MASTER.m

Note: This code requires the Statistics Toolbox
