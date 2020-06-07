res←wa angularAverage q;esta;est;p;err;tp;d;selected;w;tau;exp;log

⍝ Some usage examples:
⍝ weights ← 1÷3 × 3⍴1 		⍝ generate a vector of weights equal to 0.(3)
⍝ angles ← (?3⍴0)×○2 		⍝ generate angle samples in the 0..2pi interval
⍝ angles ← d2r (0 15 90)		⍝ store the 0 15 90 degree angles in a radian vector
⍝ wa←weights ,[0.5] angles 	⍝ store some sample weight-angle pairs
⍝ call  wa angularAverage 1 to compute the L1 weighted average/median OR
⍝ call  wa angularAverage 2 to compute the L2 weighted average/mean

 esta←(wa[2;]+.×wa[1;])÷+/wa[1;]           ⍝ compute initial angle estimate
 esta←(○2*?0)                              ⍝ or initialize it to a random 0..2pi angle
 est←*(esta×0J1)                           ⍝ and its corresponding S1 point
 exp←{⍺×(*⍵)}                              ⍝ the exp map at point ⍺ with point argument ⍵ (on S1)
 log←{⍟((+⍺)×⍵)}                           ⍝ the log map at point ⍺ with tangent vector argument ⍵ (on S1)
 p←*(wa[2;]×0J1)                           ⍝ convert angles to complex point samples
 err←○2
 iter←0
 :While (|err)>0.000001
 :AndIf iter<100
     iter+←1
     tp←est log p                          ⍝ find the tangent vectors in the space relative to the est point
     d←|tp                                 ⍝ compute tangent vector lengths
     selected←⍸(d>0.000001)                    ⍝ due to proximity to sample points, some samples must be skipped
     w←wa[1;selected]×(d[selected]*(q-2))  ⍝ compute the weighted Lq average gradient in the tangent plane
     tau←(+/(w×tp[selected]))÷(+/w)        ⍝ store it in tau
     est←est exp tau                       ⍝ project tau via the exp map at the previous estimate onto the manifold
     est÷←|est                             ⍝ normalize the estimate to counteract numerical errors
     err←(12○est)-esta                     ⍝ compute the angular deviation as an error
     esta←esta+⎕←err                       ⍝ also debug print the error as the iterations progress
 :EndWhile
 res←(12○est)                              ⍝ return the phase of the complex number
