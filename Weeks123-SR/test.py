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

print(tp,',', xp)

empty = []
a = [[0,1], [1,2]]
for l in a:
    empty.append(l)
print(empty)

empty_new = [empty]
empty_new.append(empty)
empty_new.append([[0,1]])
print(empty_new)


a = np.array([0.3, 0.4, 0])
b = np.array(a)
print(b)
print(0.3*b)
print(np.linalg.norm(a))

print(np.linalg.norm(a)*b)

if type(a) == float or int:
    print('y')
else:
    print('n')

x = [0, 1]
print(len(x))

v = [6, 7]

v_array = np.array(v) if not isinstance(v, (float, int)) else np.array([v])

print(v_array)