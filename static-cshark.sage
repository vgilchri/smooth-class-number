# we code our variant of the static DH attack for use on CSI-SharK

load("csidh.sage")

proof.all(False)

# need action function not for exponents but for taking a generator of the class group to an exponent
# to do so: find generator of class group (maybe will be expo = [1,1,...1])
# then to call g^x simply call action(base, [x,x,...x])

ls = [3,5,7]
p = 4 * prod(ls) - 1
N = 27
c = 3
v = 9
#g = 2

def csi_action(pub,exp):
	# here we assume \ell_1 generates the class group, but should verify later
	# we want to create the action from csi-fish that takes a class group gen to an exponent
	# we do so using the csidh group action
	return action(pub,[exp,0,0])

def choose_keys(N,c): # determine a secret key for csi-shark, and the necessary public info to run the attack
	# need size of class group N, and factor c of N
	sk = Integers(N).random_element()
	sk = Integers()(sk)
	E1 = csi_action(base,sk)
	Ej = csi_action(base, Integers()(sk*c))
	return sk,E1,Ej

def static_attack(E0, E1, Ej,v,c,N):
	#Ej is [(g^sk)^c]E0, taken from public key
	# E is [g^sk]E0, taken from public key
	# N = c * v
	# g generator of Z_N
	# first BSGS search-------------------
	# build subgroups
	m = ceil(sqrt(v))
	C = []
	for i in range(m):
		ind1 = i*c
		entry1 = csi_action(E0,ind1)
		C.append(entry1)
	V = []
	for j in range(m):
		ind2 = -1*m*j*c
		entry2 = csi_action(Ej,ind2)
		V.append(entry2)
	# look for collision
	for i in range(m):
		for j in range(m):
			if C[i] == V[j]:
				y=(i+j*m) % v
				break
	print("y is ",y)
	test = csi_action(base,y*c*2)
	actual = csi_action(base,c*sk*2)
	print("new E_jk is ",test)
	print("E_jk should be ",actual)
	# second BSGS search------------------
	# build subgroups
	n = ceil(sqrt(c))
	Q = []
	for i in range(n):
		ind1 = i*v
		entry1 = csi_action(E0,ind1)
		Q.append(entry1)
	W = []
	Ey = csi_action(E1,-1*y)
	for j in range(n):
		ind2 = -1*v*j*n
		entry2 = csi_action(Ey,ind2)
		W.append(entry2)
	# look for collision
	for i in range(n):
		for j in range(n):
			if Q[i] == W[j]:
				t=(i+j*n) % c
				break
	# putting it together----------------
	z = t*v + y 
	return z % (N)


#count class group number : brute force
# for p = 419, N = 27 = 3^3
#counter = 1
#end = action(base,[1,0,0])
#while base != end:
#	end = action(end,[1,0,0])
#	counter += 1
#print(counter)


sk,E1,Ej = choose_keys(N,c)

sk1 = static_attack(base, E1, Ej,v,c,N)

print(sk)

print(sk1)









