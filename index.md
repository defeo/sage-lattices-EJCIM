---
layout: index
title: Linear algebra and lattice reduction in Sage
plugins: ['mathjax', 'highlight']
---

# Linear algebra and lattice reduction in Sage 

Tutorial for the EJCIM 2014 <https://ejcim2014.greyc.fr/>.

## Linear algebra

Sage provides native support for working with matrices over any
commutative or non-commutative ring. The parent object for a matrix is
a matrix space ``MatrixSpace(R, n, m)`` of all $$n\times m$$ matrices
over a ring $$R$$.

To create a matrix, either use the ``matrix(...)`` function

~~~
sage: matrix([[1,2],[3,4]])
[1 2]
[3 4]
~~~

or create a matrix space using the ``MatrixSpace`` command and coerce
an object into it

~~~
sage: M = matrix([[1,2],[3,4]])
sage: M
[1 2]
[3 4]
~~~

Read about more ways of creating matrices here:
<http://www.sagemath.org/doc/reference/matrices/sage/matrix/constructor.html>
and here
<http://www.sagemath.org/doc/reference/matrices/sage/matrix/matrix_space.html>.

Matrices also act on row vectors, which you create using the
``vector(...)`` command or by making a
[``VectorSpace``](http://www.sagemath.org/doc/reference/modules/sage/modules/free_module.html)
and coercing lists into it.

~~~
sage: M * vector([1,1])
(3, 7)
~~~

The natural action of matrices on row vectors is from
the right. Sage currently does not have a column vector class (on
which matrices would act from the left), but this is planned.


Sage has quite flexible ways of extracting elements or submatrices
from a matrix:

~~~
sage: M[0]
(1, 2)
sage: M[0,1]
2
sage: M[:,1]
[2]
[4]
sage: _.parent()
Full MatrixSpace of 2 by 1 dense matrices over Integer Ring
~~~

You can read more on indexing here :
<http://www.sagemath.org/doc/reference/matrices/sage/matrix/docs.html#indexing>

And, obviously, you can compute the LLL and BKZ reduced forms of matrices:

~~~
sage: M.LLL()
[1 0]
[0 2]
sage: M.BKZ()
[1 0]
[0 2]
~~~

The full documentation on matrices is here :
<http://www.sagemath.org/doc/reference/matrices/index.html>. Be sure
to have read the official
[tutorial on linear algebra](http://www.sagemath.org/doc/tutorial/tour_linalg.html)
before going any further.

The [Sagebook](http://sagebook.gforge.inria.fr/) also has various
chapters on linear algebra, in particular §2.4 and §10.

#### Exercise

1. Find a basis of the solution space of the homogenous linear system
   
   $$
   A = \begin{pmatrix}
   2 & -2 & -2 & 2 & 0\\
   0 & 1/2 & -1/2 & 1 & 0\\
   0 & 0 & 1 & -1 & 0\\
   2 & -2 & -1 & 1 & -1
   \end{pmatrix}.
   $$

2. Find a basis of the space spanned by the columns of $$A$$.

## Solving Knapsacks

Before moving to lattice reduction, we study some very special
instances of lattices: those arising from *knapsack problems*.  This
section is inspired by this
[tutorial on linear programming](http://combinat.sagemath.org/doc/thematic_tutorials/linear_programming.html#linear-programming).

The *Knapsack* problem is the following: given a collection of items
having both a weight and a *usefulness*, we would like to fill a bag
whose capacity is constrained while maximizing the usefulness of the
items contained in the bag (we will consider the sum of the items'
usefulness). For the purpose of this tutorial, we set the restriction
that the bag can only carry a certain total weight.

To achieve this, we have to associate to each object $$o$$ of our
collection $$C$$ a binary variable ``taken[o]``, set to 1 when the
object is in the bag, and to 0 otherwise. Therefore, we are trying to
solve the following Mixed Integer Linear Program (MILP)

$$
\text{Max: }  \sum_{o \in L} \text{usefulness}_o \times \text{taken}_o\\
\text{Such that: }  \sum_{o \in L} \text{weight}_o \times \text{taken}_o \leq C
$$

We will restrict to the case *usefulness = weight*, which is also
known as the *subset sum* problem. Using Sage, we will give to our
items a random weight:

    sage: C = 3659
    sage: L = ["pan", "book", "knife", "gourd", "flashlight"]
    sage: L.extend(["random_stuff_" + str(i) for i in range(10)])
    sage: set_random_seed(685474)
    sage: weights = [randint(100,1000) for o in L]

We can now define the MILP itself

    sage: p = MixedIntegerLinearProgram(maximization=True)
    sage: taken = p.new_variable(binary=True)
    sage: p.add_constraint(p.sum(w * taken[o] for o,w in zip(L, weights)), max=C)
    sage: p.set_objective(p.sum(w * taken[o] for o,w in zip(L, weights)))
    sage: p.solve()
	3659.0
    sage: taken = p.get_values(taken)

The solution found is (of course) admissible

    sage: sum(w * taken[o] for o,w in zip(L, weights))
    3659.0

Should we take a flashlight?

    sage: taken["flashlight"]
    1.0

Wise advice. Based on purely random considerations.

The module
[`numerical.knapsack`](http://combinat.sagemath.org/doc/reference/numerical/sage/numerical/knapsack.html)
gives a nice wrapper for MILPs used to solve knapsack instances. The
same problem could have been solved by

	sage: from sage.numerical.knapsack import knapsack
	sage: knapsack(zip(costs,costs,L), max=C)
	[3659.0,
	 [(952, 952, 'flashlight'),
	  (767, 767, 'random_stuff_1'),
	  (286, 286, 'random_stuff_2'),
	  (152, 152, 'random_stuff_3'),
	  (806, 806, 'random_stuff_4'),
	  (696, 696, 'random_stuff_6')]]

Knapsack problems can be restated in terms of finding short vectors in
a lattice. For the sake of simplicity, let's suppose that we know
there is a solution to the problem that reaches the objective cost
$$C$$ (this is the case, for example, in the cryptanalysis of the
Merkle-Hellman system). Given items with weights $$w_0,\dots,w_n$$, we
want to find binary variables $$v_0,\dots,v_n$$ such that

$$\sum_i v_i w_i = C.$$

This can be restated as a lattice problem: consider the lattice

$$Λ = \begin{pmatrix}
1 & 0 & \cdots & 0 & w_0\\
0 & 1 & \cdots & 0 & w_1\\
\vdots &&&& \vdots\\
0 & 0 & \cdots & 1 & w_n\\
0 & 0 & \cdots & 0 & -C
\end{pmatrix},$$

then $$a=(v_0,\dots,v_n,0)$$ is a short vector in $$Λ$$. If $$a$$ is
the shortest vector, we may be lucky enough to find it using a lattice
reduction algorithm, such as LLL.


#### Exercice

1. Using the constructor `block_matrix`, construct the matrix
   $$M$$ associated to the knapsack problem above.

2. Using `LLL`, find a short binary vector in $$Λ$$. Verify that it is
   a solution to the knapsack problem.


## SIS and lattice reduction

The *Short Integer Solution* (SIS) problem is the following. Given
parameters $$n$$, $$m>2n$$ and $$q≈\mathrm{poly(n)}$$, and an $$n×m$$
matrix $$A$$ with entries uniformly distributed in $$ℤ/qℤ$$, find a
vector $$s\inℤ$$ of small norm such that

$$As=0 \mod q.$$

The SIS problem is closely related to LWE, and it is the basis of many
cryptosystems. SIS instances can be solved via lattice
reduction. Define the *dual lattice* of $$A$$

$$
Λ^\bot(A) = \{y ∈ ℤ^m \mid Ay = 0 \bmod q\},
$$

then $$s$$ is a short vector in $$Λ^\bot(A)$$.


#### Exercise

For simplicity, we are going to work with $$q$$ prime.

1. Create a SIS matrix with parameters $$10,20,1009$$, using the
   following code
   
   ~~~
   sage: n, m, q = 10, 20, 1009
   sage: set_random_seed(685474)
   sage: A = random_matrix(Zmod(q),10,20)
   ~~~
   
   Observe that $$A$$ has coefficients in $$ℤ/qℤ$$ (you can check by
   printing `A.parent()`). We are going to need matrices with
   coefficients in $$ℤ$$ in order to apply lattice reduction
   algorithms. At any time, you can obtain a copy of `A` with a
   different base ring by using `A.change_ring()`, and a copy of `A`
   with coefficients in $$ℤ$$ by using `A.lift()`.
   
2. Compute a basis of the null space of `A` in $$(ℤ/qℤ)^m$$.

3. Use this matrix to define a basis of $$Λ^\bot(A)$$ (be careful:
   this is a lattice in $$ℤ^m$$, not in $$(ℤ/qℤ)^m$).

4. Compute an LLL-reduced basis of $$Λ^\bot(A)$$. Print the norms of
   the basis vectors.

5. You have noticed that there is no vector significantly shorter in
   $$Λ^\bot(A)$$. This is no surprise, as $$A$$ is a random
   matrix. So, how does one construct an SIS instance that looks like
   a random matrix?
   
   The trick is very simple: choose a vector $$x$$ of small norm
   (e.g., a uniformly distributed binary vector), and then compute
   
   $$A' = (A \mid -Ax).$$
   
   Then $$A'$$ is an SIS instance, and we can prove that its entries
   are *almost* uniformly distributed (this is a consequence of the
   famous leftover hash lemma).
   
   Using this technique, construct an SIS instance.

6. Compute an LLL-reduced basis of $$Λ^\bot(A')$$. Can you find the
   short vector?


## LWE

We finally come to LWE. The reductions from LWE to lattice problems
are more technical than the previous ones, so we'll better leave them
off. We will instead have a look at what is included in the Sage
library concerning it.

Sage has some crypto modules built in, mainly for educational
purposes. They are listed here:
<http://www.sagemath.org/doc/reference/cryptography/index.html>.  Of
special interest to us are the
[module on lattice instance generation](http://www.sagemath.org/doc/reference/cryptography/sage/crypto/lattice.html)
and the one on
[LWE oracles](http://www.sagemath.org/doc/reference/cryptography/sage/crypto/lwe.html).

1. Read the documentation of the two modules.

2. Generate some random lattice instances and compare the results of
   LLL and BKZ reduction. How far can you go before the algorithms
   start eating all resources up?

3. Using the `DiscreteGaussianSampler`, implement the Regev and
   dual-Regev cryptosystems.
