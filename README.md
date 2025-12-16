# Partial Predictive Information Dynamics for Music 

Partial Predictive Information Dynamics for Music (ppidyom) is a new attempt at recreating (and expanding upon) [IDyOM](https://www.marcus-pearce.com/idyom/),
similar to the [ppm](https://github.com/pmcharrison/ppm) project.

---

Partial predictive modeling (PPM) has been applied widely to musical data using the IDyOM model.
However, this model tends to be a bit hard to implement, and can be inscrutable and opaque (at least to me).
Our goal is to make PPM for music faster and easier in R.
The [ppm](https://github.com/pmcharrison/ppm) is a similar idea, but our approach is more flexible and much faster.



## Algorithm


PPidyom works by doing lots of counting.
Consider this sequence:

```
I IV V I I IV V vi I ii IV V I I IV V I
```


We can chop this up into N-grams of various lengths---here we'll just show 1 through 3:


#### 1 gram

```
I IV V I I IV V vi I ii IV V I I IV V I
```

#### 2 gram

```
I IV   
  IV V 
     V I
       I I
         I IV
           IV V
              V vi
                vi I
                   I ii
                     ii IV
                        IV V
                           V I
                             I I
                               I IV
                                 IV V
                                    V I
```

#### 3 gram

```
I IV V  
  IV V I  
     V I I 
       I I IV 
         I IV V
           IV V vi
              V vi I
                vi I ii
                   I ii IV
                     ii IV V
                        IV V I
                           V I I
                             I I IV
                               I IV V
                                 IV V I
```

Now, for each length of N-gram, we'll walk through the sequence and count how many times each gram has been observed *before*.
(We'll count at the last token of each gram.) 

#### 1 gram

```
I IV V I I IV V vi I ii IV V I I IV V I
0 0  0 1 2 1  1 0  3 0  2  2 4 5 3  3 6
```

By the time we get to the end, `I` has been observed seven times (it says size because it had been observed six times *before* that final one).

#### 2 gram

```
  
I IV 0
  IV V 0
     V I 0
       I I 0
         I IV 1
           IV V 1
              V vi 0
                vi I 0
                   I ii 0
                     ii IV 0
                        IV V 2
                           V I 1
                             I I 1
                               I IV 2
                                 IV V 3
                                    V I 2
```

By the time we get to the end, we've observed `IV V` the most (four times).
With the final `V I` cadence, we have heard that progression two times previously.

#### 3 gram

```
I IV V 0   
  IV V I 0  
     V I I 0 
       I I IV 0
         I IV V 1
           IV V vi 0 
              V vi I 0
                vi I ii 0
                   I ii IV 0
                     ii IV V 0
                        IV V I 1
                           V I I 1
                             I I IV 1
                               I IV V 2
                                 IV V I 2
```


### Getting probabilities

Ok, let's line up the counts we've got so far:

```
        I IV V I I IV V vi I ii IV V  I  I  IV V  I
0-gram: 0 1  2 3 4 5  6 7  8 9  10 11 12 13 14 15 16
1-gram: 0 0  0 1 2 1  1 0  3 0  2  2  4  5  3  3  6
2-gram:   0  0 0 0 1  1 0  0 0  0  2  1  1  2  3  2
3-gram:      0 0 0 0  1 0  0 0  0  0  1  1  1  2  2
```

We'll call this our "count matrix."
Given our count matrix, what can we calculated?

Let's think about it like this: as we listen through the piece, one chord a time, we can remember how many times each N-gram has occurred before (that's what our counts are!).
When we get to penultimate `V` chord, we could stop and ask, based on what we've seen how probable is it that the next (final) chord will be `I`?
Well, at this point we've heard `V` three times before.
In the past, we've heard `V I` twice, and `V vi` once.
So, based on our past experience, there would be a 2/3 probability of `V I` and a 1/3 probability of `V vi`.
Thus, we can get the 1st-order dynamic conditional probabilities of our sequence quite easily:

$$
P(chord_i \| chord_{i-1}) = \frac{Count2Gram(chord_i)}{Count1Gram(chord_{i-1})}\\
\intertext{where}
i = index
$$


We can apply the exact same formula with the three-grams and two-grams to get the 2nd-order conditional probability:

$$
P(chord_i \| chord_{i-1, i-2}) = \frac{Count3Gram(chord_i)}{Count2Gram(chord_{i-1})}\\
\intertext{where}
i = index
$$

For example, for the final `I`, the probability of of that chord being `I` given the previous two chords is also 2/3.
Why? `IV V I` has been heard twice before, while `IV V` has been heard three times before.


> Notice that if you want to know the posterior probability of a chord, taking into account that you've now seen it, just add 1 to the numerator and denominator.

### Why use the "count matrix"?

The nice things about the count-matrix approach include that 

1) Have more data? It is easy to just add more counts! This allows us to flexibility model and update prior knowledge.
2) Similarly, we can easily add "escape probabilities" and other things we need to the table.
3) If we've never seen a N-gram before, we can easily "back off" to the smaller N-grams.
4) It is fast. Using the `data.table` package, we can compute the count matrix very fast, even for large sequences.


