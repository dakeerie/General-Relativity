import numpy as np
import sympy

def function(x):
    print(x)

matrix = np.array([[1,2,3], [1,2,3]])

vector = np.array([[1,2,3]])
column = vector.transpose()

print(matrix, vector)

print(np.shape(vector))
print(np.shape(column))


print(np.matmul(vector, column))

print(np.identity(3))

sq = np.array([[1, 2, 0], [2, 1, 0], [0, 0, 2]])
print(np.linalg.inv(sq))
print(np.matmul(sq, np.linalg.inv(sq)))
print(np.matmul(sq, np.identity(3), np.linalg.inv(sq)))


print(sq @ np.linalg.inv(sq) @ sq)

x = ('apple', 'banana')
for index, fruit in enumerate(x):
    print(index)


x, y = sympy.symbols('x y')
gamma, v = sympy.symbols('gamma v')
lt = np.array([[gamma, -gamma*v], [-gamma*v, gamma]])

tp, xp = lt @ np.array([x,y])

print(tp, xp)