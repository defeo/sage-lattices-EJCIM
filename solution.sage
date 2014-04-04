## Knapsack

C = 3659
L = ["pan", "book", "knife", "gourd", "flashlight"]
L.extend(["random_stuff_" + str(i) for i in range(10)])
set_random_seed(685474)
weights = [randint(100,1000) for o in L]

#### we build the knapsack lattice
M = block_matrix(2,2,[identity_matrix(len(weights)), 0, matrix(weights), -C]).transpose()
#### and find a solution vector in the last row of the matrix
sol = M.LLL()[-1]
#### we verify that it is indeed a solution
assert all(map(lambda x: x in (0, 1), sol))
assert sol.inner_product(vector(weights + [0])) == C


## SIS

n, m, q = 10, 20, 1009
set_random_seed(685474)
A = random_matrix(Zmod(q),10,20)

#### The null space
Anull = A.right_kernel_matrix()
#### The dual lattice
LambdaA = block_matrix(2, 1, [Anull.lift(), q*identity_matrix(20)])
reduced = LambdaA.LLL()
print min(filter(lambda x:x, [r.norm() for r in reduced]))

#### An SIS instance
x = vector(ZZ, [randint(0,1) for i in range(20)])
A1 = block_matrix(1, 2, [A,matrix(-A*x).transpose()])
#### Its null space
A1null = A1.right_kernel_matrix()
#### Its dual lattice
LambdaA1 = block_matrix(2, 1, [A1null.lift(), q*identity_matrix(21)])
reduced1 = LambdaA1.LLL()
norms = [r.norm() for r in reduced1]
print min(filter(lambda x:x, norms))
assert 2*sqrt(2) in norms
assert reduced1[norms.index(2*sqrt(2))][:-1] == x
