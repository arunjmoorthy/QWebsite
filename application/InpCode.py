import numpy.linalg
import math
import random

##sets the dimension throughout
##always is a prime
q = 3


##converts integer n to a q-ary value with pts digits
def inttoarr(n, pts):
    vals = numpy.zeros(pts)
    for k in range(pts):
        vals[k] = int(n / math.pow(q, pts - 1 - k))
        n = n - vals[k] * math.pow(q, pts - 1 - k)
    return vals
## q=2
##int n - 8 and points = 4
##list of q-ary -  1000
## q=3
## n - 8 points = 4
## 2(3) + 2(1) -- 0022

##returns the full symp group
##takes a stabilizer "generator" set c and prints the full group
def sympgroup2(c):
    sympgroup = []
    powers = numpy.zeros(len(c))
    for i in range(int(math.pow(q, len(c)))):
        powers = inttoarr(i, len(c))
        out = numpy.zeros(len(c[0]))
        for j in range(len(c)):
            temp = numpy.zeros(len(c[0]))
            for k in range(len(c[0])):
                temp[k] = powers[j] * c[j][k]
            out = numpy.add(out, temp)
            out = out % q
        sympgroup.append(out.tolist())
    return sympgroup

##takes generators and generators stabilizer group
##takes each generator - composes it w other generators
##adds copies of the generators together - forms full stabilizer group
## q=3
## xx, zz
## 2 generators - 0,1,2
## if xx = 1, zz = 0 --> get xx
## compose xx with zz --> x*z tensor x*z
## xx ^ 0 means Identity tensor identity
## xx ^ 1 means one composition of xx - zz ^ 2 --> x*z^2 tensor x*z^2
## purpose of generating full stabilizer group --> need to know if error is in stabilizer group
## errors are 2n length vectors
## commutator btw error and all generators - if error is in group, information protected
## output is k values - commutator between error and all gen
## output vector - all in that vector are congruent to 0 then its undetectable (means stabilizer code doesn't recognize error)
## that error could be in stabilizer code
## if error doesn't commute w all stabilizer generators, if it is a low weight error, you can know what the error was (and undo it)
## can decompose error into generalized Pauli
## 2n vectors can consist of all errors possible

## have information - doing sumthn with it
## environmental noise gets in and messes up information a bit
## physical error is some unitary - can go to generalized Pauli's
## w generalized Pauli's, you can map to 2n vector w symplectic representation
## ex. single Pauli error --> commutator between code and error
## syndrome values correspond to a Pauli, so by looking at them u can know what happened (what the error was)
## which syndromes correspond to which errors


##not updated to qudits
##converts symp to pauli
## works only w qubits
def symptop(symp):
    p = ''
    n = int(len(symp) / 2)
    for i in range(n):
        if (symp[i] == 1 and symp[i + n] == 1):
            p = p + 'Y'
        elif (symp[i] == 1):
            p = p + 'X'
        elif (symp[i + n] == 1):
            p = p + 'Z'
        else:
            p = p + 'I'
    return p



##not updated to qudits
##converts pauli to symp
def ptosymp(pstr):
    symp = numpy.zeros(2 * len(pstr))
    for i in range(len(pstr)):
        if (pstr[i] == 'X'):
            symp[i] = 1
        if (pstr[i] == 'Z'):
            symp[i + len(pstr)] = 1
        if (pstr[i] == 'Y'):
            symp[i] = 1
            symp[i + len(pstr)] = 1
    symp = symp.astype(int)
    return symp.tolist()


##verifies that the elements of c satisfy the requirements for stabilizer code
def verify(c):
    # returning 0 means it's valid [all commute], anything else is invalid
    valid = 0
    for i in range(len(c)):
        for j in range(len(c)):
            valid = valid + comm(c[i], c[j])
    return valid

## assumes that there is a subset of stabilizer code (not just generators)
## Checks whether all pairs of generators commute w each other


##computes the commutator of two symp paulis
## can be btw error and generator or just generators
## with mod q
def comm(p1, p2):
    if (not (len(p1) == len(p2))):
        print("dim error")
        return
    comm = 0
    n = int(len(p1) / 2)
    # print(n)
    for i in range(n):
        # print(p1[i])
        # print(p2[(i+len(p1))%len(p1)])
        comm = comm + (p1[i] * p2[i + n] - p1[i + n] * p2[i])
    return comm % q


##computes the commutator of two symp paulis
## w o mod q
def comminf(p1, p2):
    if (not (len(p1) == len(p2))):
        print("dim error")
        return
    comm = 0
    n = int(len(p1) / 2)
    # print(n)
    for i in range(n):
        # print(p1[i])
        # print(p2[(i+len(p1))%len(p1)])
        comm = comm + (p1[i] * p2[i + n] - p1[i + n] * p2[i])
    return comm


##this will go through the group and check if any sum is repeated
##takes a set of Pauli and removes any redundant ones
##linear combination of another element in the matrix
def makegens2(c):
    powers = numpy.zeros(len(c))
    newc = c.copy()
    removed = numpy.zeros(len(c))
    for i in range(int(math.pow(q, len(c)))):
        powers = inttoarr(i, len(c))
        out = numpy.zeros(len(c[0]))
        tally = 0
        for l in range(len(c)):
            tally = tally + ((powers[l] > 0) and removed[l])
        if (tally > 0):
            continue
        for j in range(len(c)):
            temp = numpy.zeros(len(c[0]))
            for k in range(len(c[0])):
                temp[k] = powers[j] * c[j][k]
            out = numpy.add(out, temp)
        out = out % q
        out = out.astype(int)
        out = out.tolist()
        if ((out in newc) and not (math.log(i, q) == int(math.log(i, q)))):
            newc.remove(out)
            removed[c.index(out)] = 1
    return newc


##computes the pauli weight of a symp pauli
## add 1 if i and n+i are either non-zero - add 0 if both are 0
def pweight(s):
    weight = 0
    n = int(len(s) / 2)
    for k in range(n):
        weight = weight + ((s[k] > 0) or (s[k + n] > 0))
    return weight


##returns the distance of the code through brute-force
##assumes gens are gens
##assumes nondegenerate code
##how to build up by Pauli weight? --use aux function
def dist(c):
    dmin = len(c[0])
    symp = sympgroup2(c)
    for i in range(int(math.pow(q, len(c[0])))):
        if (math.log(i, q) == int(math.log(i, q))):
            print(i)
        err = inttoarr(i, len(c[0]))
        err = err.tolist()
        temp = 0
        for j in range(len(c)):
            temp = temp + (comm(err, c[j]) % q)
        ##if it has a nonzero syndrome, skip it
        if (not (temp == 0)):
            continue
        tem = (err in symp)
        if (not (tem)):
            dmin = min(dmin, pweight(err))
    return dmin


def dist2(c):
    dmin = len(c[0])
    symp = sympgroup2(c)
    found = 0
    weight = 1
    while (not (found)):
        errs = hammstrs(int(len(c[0]) / 2), weight)
        if (errs == []):
            found = 1
            print("distance is larger than 3")
            return
        for i in range(len(errs)):
            temp = 0
            for j in range(len(c)):
                temp = temp + (comm(errs[i], c[j]) % q)
            if (not (temp == 0)):
                continue
            tem = (errs[i] in symp)
            if (not (tem)):
                dmin = weight
                found = 1
                return dmin
        weight = weight + 1
    return dmin


##generates all Hamming strings of length n with weight w
## Hamming weight is number of non-zero entries
## input is number of registers
## outputs all binary strings with Hamming weight w
def hammstrs(n, w):
    # we will set t=q later
    t = q
    allstrs = []
    if (w == 0):
        temp = []
        for j in range(2 * n):
            temp.append(0)
        return temp
    # total shift values
    if (w == 1):
        print("checking w=1")
        for j in range(n):
            for q1 in range(t):
                for q2 in range(t):
                    temp = []
                    for k in range(j):
                        temp.append(0)
                    temp.append(q1)
                    for k in range(n - j - 1):
                        temp.append(0)
                    for k in range(j):
                        temp.append(0)
                    temp.append(q2)
                    for k in range(n - j - 1):
                        temp.append(0)
                    # print(temp)
                    allstrs.append(temp)
        return allstrs
    if (w == 2):
        ##loops need some fixing, but patched for now
        for j1 in range(n - 1):
            for j2 in range(j1, n - 1):
                for q1 in range(t):
                    for q2 in range(t):
                        for q3 in range(t):
                            for q4 in range(t):
                                temp = []
                                for k in range(j1):
                                    temp.append(0)
                                temp.append(q1)
                                for k in range(j2 - j1 - 1):
                                    temp.append(0)
                                temp.append(q3)
                                for k in range(n - j2 - 2):
                                    temp.append(0)
                                for k in range(j1):
                                    temp.append(0)
                                temp.append(q2)
                                for k in range(j2 - j1 - 1):
                                    temp.append(0)
                                temp.append(q4)
                                for k in range(n - j2 - 2):
                                    temp.append(0)
                                if (len(temp) == 2 * n):
                                    allstrs.append(temp)
        return allstrs
    if (w == 3):
        for j1 in range(n - 2):
            for j2 in range(j1, n - 2):
                for j3 in range(j2, n - 2):
                    for q1 in range(t):
                        for q2 in range(t):
                            for q3 in range(t):
                                for q4 in range(t):
                                    for q5 in range(t):
                                        for q6 in range(t):
                                            temp = []
                                            for k in range(j1):
                                                temp.append(0)
                                            temp.append(q1)
                                            for k in range(j2 - j1 - 1):
                                                temp.append(0)
                                            temp.append(q3)
                                            for k in range(j3 - j2 - 1):
                                                temp.append(0)
                                            temp.append(q5)
                                            for k in range(n - j3 - 3):
                                                temp.append(0)
                                            for k in range(j1):
                                                temp.append(0)
                                            temp.append(q2)
                                            for k in range(j2 - j1 - 1):
                                                temp.append(0)
                                            temp.append(q4)
                                            for k in range(j3 - j2 - 1):
                                                temp.append(0)
                                            temp.append(q6)
                                            for k in range(n - j3 - 3):
                                                temp.append(0)
                                            if (len(temp) == 2 * n):
                                                allstrs.append(temp)
        return allstrs
        # continue in this way up to say w=5
    # this goes by Pauli weight
    # recursion is the way to go...
    return []


##returns the multiplicative inverse of a number mod q
def minv(x, q):
    if (x == 0):
        return -1
    for i in range(1, q):
        if (int((i * x) % q) == 1):
            return i
    return "error"

##this performs the effect of conjugation by dft on column col
def dft(c, col):
    temp = 0
    for i in range(len(c)):
        temp = c[i][col]
        c[i][col] = c[i][col + int(len(c[0]) / 2)]
        c[i][col + int(len(c[0]) / 2)] = (-1 * temp) % q
    return c

## performs a discrete fourier transform
## switch the roles of x and z operators on physical qudit
## tool to put a code into canonical form


##this performs a swapping of qudits reg1 and reg2
def regswap(c, reg1, reg2):
    return


##finite field reduced row echelon form--used for putting the code into canonical form
##should be good, sans oversights and mentioned errors
def ffrref(c):
    outc = numpy.zeros([len(c), len(c[0])])
    for i in range(len(c)):
        x = c[i][i]
        inv = minv(x, q)
        # print(inv)
        if (inv == -1):
            # print(c)
            # print("severe issue")
            temp = outc.astype(int)
            return temp.tolist()
        ##        if(inv==-1):
        ##            ##need to use DFT to fix this
        ##            c=dft(c,i)
        ##            x=c[i][i]
        ##            inv=minv(x,q)
        ##            ##next we need to do register swaps. find row>i with either c[row][i]>0 or c[row][i+n]>0
        ##            if(inv==-1):
        ##                return "crash due to support in wrong place"
        for j in range(len(c[0])):
            outc[i][j] = (inv * c[i][j]) % q
        for k in range(len(c)):
            if (not (k == i)):
                for j in range(len(c[0])):
                    outc[k][j] = (c[k][j] - c[k][i] * c[i][j]) % q
    outc = outc.astype(int)
    return outc.tolist()


##still not fully working
####do we need to make an empty numpy array again?
def symptocon(c):
    n = len(c[0]) // 2
    for i in range(len(c)):
        if ((c[i][i] == 0)):
            stop = 0
            for j in range(i, len(c)):
                if (stop > 0):
                    continue
                f = n + i
                if (not (c[j][i] == 0)):
                    temp = c[i]
                    c[i] = c[j]
                    c[j] = temp
                    stop = 1
                # need to do the Fourier Transform
                elif (not (c[j][f] == 0)):
                    temp = c[i]
                    c[i] = c[j]
                    c[j] = temp
                    c = dft(c, i)
                    stop = 1
    return c


##generates a random stabilizer code over n qudits with k generators
def genrandoms(n, k):
    currcode = []
    toadd = []
    #creates a random pauli
    for m in range(2 * n):
        toadd.append(random.randint(0, q - 1))

    currcode.append(toadd)
    counter = 0
    # makes sure the code doesn't get stuck in a bad seed forever
    #if the code doesn't have enough generators yet, then keep trying
    while (len(currcode) < k and counter < 50000):
        toadd = []
        counter = counter + 1
        # if(counter%100==0):
        #    print(counter)

        #generates another random Pauli
        for m in range(2 * n):
            toadd.append(random.randint(0, q - 1))
        include = 1
        #assuming good to add to list
        #then checks if it is fine
        for i in range(len(currcode)):
            if (not (comm(currcode[i], toadd) % q == 0)):
                include = 0
        #if it is good to go, include it in the currcode
        if (include):
            temp = currcode.copy()
            temp.append(toadd)
            if (len(makegens2(temp)) > len(currcode)):
                currcode.append(toadd)
    return currcode


def verifyinf(c):
    # returning 0 means it's valid [all commute], anything else is invalid
    valid = 0
    for i in range(len(c)):
        for j in range(len(c)):
            valid = valid + comminf(c[i], c[j])
    return valid


def lmatrix(c):
    #    l = [[0] * len(c)] * len(c)
    l = numpy.zeros([len(c), len(c)])
    for i in range(len(c)):
        for j in range(len(c)):
            if (i >= j):
                n = comminf(c[i], c[j])
                l[i][j] = n
    l = l.astype(int)
    l = l.tolist()
    return l


##finish c+[0_{kxk} 0_{kx(n-k)} | L_{kxk} 0_{kx(n-k)}]
##c[i][j+n]=L[i][j] 0<i,j<k
def invariant(c, l):
    n = len(c[0]) // 2
    k = len(c)
    for i in range(k):
        for j in range(k):
            c[i][j + n] += l[i][j]
    return c


def alltogetherModP(c, q, p):
    dq = dist2(c)
    c = symptocon(c)
    c = ffrref(c)
    c = invariant(c, lmatrix(c))
    dp = dist2(c)
    distances = [dq, dp]
    r = [c, distances, ""]
    if(dp >= dq):
        r[2] = "Yes"
    else:
        r[2] = "No"

    return r

def alltogetherCan(c, q, p):
    dq = dist2(c)
    c = symptocon(c)
    c = ffrref(c)
    r = [c, "", ""]

    return r

def alltogetherInv(c, q, p):
    dq = dist2(c)
    c = symptocon(c)
    c = ffrref(c)
    c = invariant(c, lmatrix(c))
    dp = dist2(c)
    distances = [dq, dp]
    r = [c, distances, ""]
    if(dp >= dq):
        r[2] = "Yes"
    else:
        r[2] = "No"

    return r

