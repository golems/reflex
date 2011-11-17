#!/usr/bin/octave

A = [0,-1;1,-.5]'

B = [0;1]

C = [1,0]

V = [1 0; 0 1]
W = 1


# kalman-bucy

ricB =  C'*inv(W)*C

P = are(A',ricB, V)

K = P * C' * inv(W)
