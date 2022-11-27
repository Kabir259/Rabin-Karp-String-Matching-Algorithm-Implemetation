import random
import math
import string

alphabet = list(string.ascii_uppercase)
d = 26


def hashbasic(a):
    lst = [('?', 0)]
    for i in range(len(a)):
        lst.append((a[i], i))
    return dict(lst)


a = hashbasic(alphabet)
# basic mapping, '?' is mapped to 0

# print(a)


# To generate random prime less than N
def randPrime(N):
    primes = []
    for q in range(2, N + 1):
        if isPrime(q):
            primes.append(q)
    return primes[random.randint(0, len(primes) - 1)]


# To check if a number is prime
def isPrime(q):
    if q > 1:
        for i in range(2, int(math.sqrt(q)) + 1):
            if q % i == 0:
                return False
        return True
    else:
        return False


# pattern matching
def randPatternMatch(eps, p, x):
    N = findN(eps, len(p))
    q = randPrime(N)
    return modPatternMatch(q, p, x)


# pattern matching with wildcard
def randPatternMatchWildcard(eps, p, x):
    N = findN(eps, len(p))
    q = randPrime(N)
    return modPatternMatchWildcard(q, p, x)


# return appropriate N that satisfies the error bounds
def findN(eps, m):
    p = m / eps
    return math.floor(2 * p * math.log2(p)) + 1


'''Proof for the claim of findN
Given - eps is the upper bound on the probability of the algorithm giving a False Positive when hash values aren't equal. 
Lets call hash of text and hash of pattern as x and y resp. (same length 'm')
False positive is reported when |x-y|%q = 0
If D = |x-y| then D is also of the same length m.
WKT x,y,D < 26^m
Any D can be factorised into multiple primes and with each prime being greater than 2, if there are 'k' no. of prime factors, we can say D >= 2^k
So k <= logD.  Hence property 1 given in the assignment is proved.

Now as D < 26^m, we can say k < mlog26 (if this was in binary then k<m).

The probability that the randomly chosen prime is one in the set specified by us {1,2...N} is 
m/pi(N) < eps
pi(N)>m/eps, let m/eps be called 'a'

CLAIM: Take N = ceiling(2a.log(a))
With this N, pi(N) >= a or m/eps or m (if eps >= 1)

Pr[False positive is reported] <= m/pi(N) <= m/eps*m = (m/m)*(1/eps) = 1/eps.
Hence Proved.

1) pi(n) > (1.26...)n/logn > n/logn (chebyshev inequality)
2) Density of primes: if we want at least k ≥ 4 primes between 1 and n, it suffices to have n ≥ 2klogk
               Proof: π(2klogk)≥ 2klogk/log(2klogk) using (1)
                      2klogk/log(2klogk)  ≥ 2klogk/[log2 + logk + log(logk) + 2] ≥ k 
                      Hence proved
                      
The no. of bits required for this operation is 2logN when N is as mentioned above
N = 2sm(log(sm)) (say, and s = 1/eps), so 2logN = 2logs + 2logm + 2log(log(sm)) + 2
2logN = O(log(s)+log(m))
s here is 1/eps
So no. of bits required are O(1/eps + logm)

Time complexity is O(1)
'''


# Return sorted list of starting indices where p matches x
def modPatternMatch(q, p, x):
    text = x
    pattern = p
    n = len(text)
    m = len(pattern)
    h = pow(d, m - 1) % q
    p = 0
    t = 0
    ans = []

    for i in range(m):
        p = (d * p + a[pattern[i]]) % q
        t = (d * t + a[text[i]]) % q

    for s in range(n - m + 1):
        if p == t:
            ans.extend([s])

        if s < n - m:
            t = (t - h * a[text[s]]) % q  # remove letter s
            t = (t * d + a[text[s + m]]) % q  # add letter s+m
            t = (t + q) % q  # make sure that t >= 0
    return ans


def modPatternMatchWildcard(q, p, x):
    text = x
    pattern = p
    n = len(text)
    m = len(pattern)
    wildcard = pattern.index('?')
    h = d ** (m - 1) % q  # this will be used as a factor which the highest digit of the hashvalue will multiplied and then subtracted from the main value so that the hash window can be rolled forward
    hPlaceholder = d ** (m - wildcard - 1) % q  # same with this as stated above but for the wildcard case
    p = 0
    t = 0
    tPlaceholder = 0
    ans = []

    for i in range(m):  # hash of pattern which skips over '?'
        if i == wildcard:
            continue
        p = (p + pow(26, m - i - 1) * a[pattern[i]]) % q

    for i in range(m):  # hash of text, ('?' is assigned 0 in the dictionary)
        t = (t + 26 ** (m - i - 1) * a[text[i]]) % q

    for i in range(m):  # modified hash of text which skips over '?' and acts as a skeleton hesh
        if i == wildcard:
            continue
        tPlaceholder = (tPlaceholder + 26 ** (m - i - 1) * a[text[i]]) % q

    if p == tPlaceholder:
        ans.extend([0])

    for s in range(1, n - m + 1):
        tPlaceholder1 = (t - h * a[text[s - 1]]) * d  # take the original hash and remove highest digit
        tPlaceholder2 = tPlaceholder1 + a[x[s + m - 1]]  # add new lowest digit
        tPlaceholder3 = tPlaceholder2 - (hPlaceholder * a[text[s + wildcard]])  # subtract from the skeleton created by the original hash the value of wildcard at appropriate location
        tPlaceholder = tPlaceholder3 % q  # assign the skeleton hash to its modded value

        t = (t - h * a[text[s - 1]]) % q  # remove letter s
        t = (t * d + a[text[s + m - 1]]) % q  # add letter s+m
        t = (t + q) % q  # make sure that t >= 0

        if p == tPlaceholder:
            ans.extend([s])

    return ans
